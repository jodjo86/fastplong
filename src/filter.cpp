#include "processor.h"
#include "seprocessor.h"

Filter::Filter(Options* opt){
    mOptions = opt;
}


Filter::~Filter(){
}

int Filter::passFilter(Read* r) {
    if(r == NULL || r->length()==0) {
        return FAIL_LENGTH;
    }

    int rlen = r->length();
    int lowQualNum = 0;
    int nBaseNum = 0;
    int gcNum = 0;
    int totalQual = 0;

    // need to recalculate read-level metrics if the corresponding filters are enabled
    if(mOptions->qualfilter.enabled || mOptions->gcContentFilter.enabled) {
        const char* seqstr = r->mSeq->c_str();
        const char* qualstr = r->mQuality->c_str();

        for(int i=0; i<rlen; i++) {
            char base = seqstr[i];
            char qual = qualstr[i];

            totalQual += qual - 33;

            if(qual < mOptions->qualfilter.qualifiedQual)
                lowQualNum ++;

            if(base == 'N')
                nBaseNum++;

            if(base == 'G' || base == 'C' || base == 'g' || base == 'c')
                gcNum++;
        }
    }

    if(mOptions->qualfilter.enabled) {
        if(lowQualNum > (mOptions->qualfilter.unqualifiedPercentLimit * rlen / 100.0) )
            return FAIL_QUALITY;
        else if(mOptions->qualfilter.avgQualReq > 0 && (totalQual / rlen)<mOptions->qualfilter.avgQualReq)
            return FAIL_QUALITY;
        else if(nBaseNum * 100 > rlen * mOptions->qualfilter.nBasePercentLimit )
            return FAIL_N_BASE;
        else if(mOptions->qualfilter.nBaseLimit !=1000000 && nBaseNum > mOptions->qualfilter.nBaseLimit)
            return FAIL_N_BASE;
    }

    if(mOptions->lengthFilter.enabled) {
        if(rlen < mOptions->lengthFilter.requiredLength)
            return FAIL_LENGTH;
        if(mOptions->lengthFilter.maxLength > 0 && rlen > mOptions->lengthFilter.maxLength)
            return FAIL_TOO_LONG;
    }

    if(mOptions->gcContentFilter.enabled) {
        double gcContent = gcNum * 100.0 / rlen;
        if(gcContent < mOptions->gcContentFilter.min || gcContent > mOptions->gcContentFilter.max)
            return FAIL_GC_CONTENT;
    }

    if(mOptions->complexityFilter.enabled) {
        if(!passLowComplexityFilter(r))
            return FAIL_COMPLEXITY;
    }

    return PASS_FILTER;
}

bool Filter::passLowComplexityFilter(Read* r) {
    int diff = 0;
    int length = r->length();
    if(length <= 1)
        return false;
    const char* data = r->mSeq->c_str();
    for(int i=0; i<length-1; i++) {
        if(data[i] != data[i+1])
            diff++;
    }
    if( (double)diff/(double)(length-1) >= mOptions->complexityFilter.threshold )
        return true;
    else
        return false;
}

vector<pair<int, int>> Filter::detectLowQualityRegions(Read* r, int windowSize, int quality) {
    vector<pair<int, int>> results;
    if(r == NULL || r->length() == 0 || windowSize <=0)
        return results;

    int l = r->length();
    const char* qualstr = r->mQuality->c_str();

    int start = 0;
    while(start + windowSize <= l) {
        int totalQual = 0;
        // preparing rolling
        for(int i=start; i<windowSize-1 && i<l; i++)
            totalQual += qualstr[i];

        int windowStart = -1;
        // find the first window with mean quality < quality
        for(int s=start; s+windowSize<l; s++) {
            if(totalQual < (33 + quality) * windowSize) {
                windowStart = s;
                break;
            }
            // roll to the new base
            totalQual += qualstr[s+windowSize];
            totalQual -= qualstr[s];
        }

        if(windowStart == -1)
            break;

        //extend the window
        int e;
        for(e=windowStart; e+windowSize<l; e++) {
            // roll to the new base
            totalQual += qualstr[e+windowSize];
            totalQual -= qualstr[e];
            if(totalQual >= (33 + quality) * windowSize) {
                break;
            }
        }
        results.push_back(make_pair(windowStart, e+windowSize-1));
        start = e + windowSize;
    }

    return results;
}

Read* Filter::keepBestReadSegment(Read* r, int windowSize, int quality, int& trimmedBases) {
    trimmedBases = 0;
    if(r == NULL || r->length() <= 0)
        return NULL;

    int readLen = r->length();
    const char* qualstr = r->mQuality->c_str();
    vector<pair<int, int>> lowQualRegions;
    if(readLen < windowSize) {
        long totalQual = 0;
        for(int i=0; i<readLen; i++)
            totalQual += qualstr[i] - 33;
        if((double)totalQual / (double)readLen < quality)
            return NULL;
        return r;
    }

    long windowQual = 0;
    for(int i=0; i<windowSize; i++)
        windowQual += qualstr[i] - 33;
    for(int start=0; start + windowSize <= readLen; start++) {
        if(start > 0) {
            windowQual -= qualstr[start - 1] - 33;
            windowQual += qualstr[start + windowSize - 1] - 33;
        }
        if((double)windowQual / (double)windowSize < quality) {
            int lowStart = start;
            int lowEnd = start + windowSize - 1;
            if(!lowQualRegions.empty() && lowStart <= lowQualRegions.back().second + 1) {
                lowQualRegions.back().second = lowEnd;
            } else {
                lowQualRegions.push_back(make_pair(lowStart, lowEnd));
            }
        }
    }
    if(lowQualRegions.empty())
        return r;

    vector<pair<int, int>> refinedLowQualRegions;
    for(size_t i=0; i<lowQualRegions.size(); i++) {
        int start = lowQualRegions[i].first;
        int end = lowQualRegions[i].second;
        while(start <= end && qualstr[start] - 33 >= quality)
            start++;
        while(end >= start && qualstr[end] - 33 >= quality)
            end--;
        if(start <= end)
            refinedLowQualRegions.push_back(make_pair(start, end));
    }
    lowQualRegions = refinedLowQualRegions;
    if(lowQualRegions.empty())
        return r;

    vector<pair<int, int>> segments;
    int lastEnd = -1;
    for(size_t i=0; i<lowQualRegions.size(); i++) {
        int start = max(0, lowQualRegions[i].first);
        int end = min(readLen - 1, lowQualRegions[i].second);
        if(start > lastEnd + 1)
            segments.push_back(make_pair(lastEnd + 1, start - 1));
        lastEnd = max(lastEnd, end);
    }
    if(lastEnd < readLen - 1)
        segments.push_back(make_pair(lastEnd + 1, readLen - 1));

    if(segments.empty())
        return NULL;

    int bestStart = -1;
    int bestEnd = -1;
    double bestMeanQual = -1.0;
    int bestLength = 0;
    for(size_t i=0; i<segments.size(); i++) {
        int start = segments[i].first;
        int end = segments[i].second;
        if(start < 0 || end >= readLen || start > end)
            continue;
        long totalQual = 0;
        for(int p=start; p<=end; p++)
            totalQual += qualstr[p] - 33;
        int len = end - start + 1;
        double meanQual = (double)totalQual / (double)len;
        if(meanQual > bestMeanQual ||
           (meanQual == bestMeanQual && len > bestLength) ||
           (meanQual == bestMeanQual && len == bestLength && start < bestStart)) {
            bestMeanQual = meanQual;
            bestLength = len;
            bestStart = start;
            bestEnd = end;
        }
    }

    if(bestStart < 0 || bestEnd < bestStart)
        return NULL;

    trimmedBases = readLen - bestLength;
    if(trimmedBases <= 0)
        return r;

    string* seq = new string(*r->mSeq, bestStart, bestLength);
    string* qual = new string(*r->mQuality, bestStart, bestLength);
    string* name = new string(*r->mName);
    name->insert(1, "best-segment-");
    string* strand = new string(*r->mStrand);
    return new Read(name, seq, strand, qual);
}

Read* Filter::trimAndCut(Read* r, int front, int tail, int& frontTrimmed) {
    frontTrimmed = 0;
    // return the same read for speed if no change needed
    if(front == 0 && tail == 0 && !mOptions->qualityCut.enabledFront && !mOptions->qualityCut.enabledTail)
        return r;


    int rlen = r->length() - front - tail ; 
    if (rlen < 0)
        return NULL;

    if(front == 0 && !mOptions->qualityCut.enabledFront && !mOptions->qualityCut.enabledTail){
        r->resize(rlen);
        return r;
    } else if(!mOptions->qualityCut.enabledFront && !mOptions->qualityCut.enabledTail){
        r->mSeq->erase(0,front);
        r->mSeq->resize(rlen);
        r->mQuality->erase(0,front);
        r->mQuality->resize(rlen);
        frontTrimmed  = front;
        return r;
    }

    // need quality cutting

    int l = r->length();
    const char* qualstr = r->mQuality->c_str();
    const char* seq = r->mSeq->c_str();
    // quality cutting forward
    if(mOptions->qualityCut.enabledFront) {
        int w = mOptions->qualityCut.windowSizeFront;
        int s = front;
        if(l - front - tail - w <= 0)
            return NULL;

        int totalQual = 0;

        // preparing rolling
        for(int i=0; i<w-1; i++)
            totalQual += qualstr[s+i];

        for(s=front; s+w<l-tail; s++) {
            totalQual += qualstr[s+w-1];
            // rolling
            if(s > front) {
                totalQual -= qualstr[s-1];
            }
            // add 33 for phred33 transforming
            if((double)totalQual / (double)w >= 33 + mOptions->qualityCut.qualityFront)
                break;
        }

        // the trimming in front is forwarded and rlen is recalculated
        if(s >0 )
            s = s+w-1;
        while(s<l && seq[s] == 'N')
            s++;
        front = s;
        rlen = l - front - tail;
    }
    // quality cutting backward
    if(mOptions->qualityCut.enabledTail) {
        int w = mOptions->qualityCut.windowSizeTail;
        if(l - front - tail - w <= 0)
            return NULL;

        int totalQual = 0;
        int t = l - tail - 1;

        // preparing rolling
        for(int i=0; i<w-1; i++)
            totalQual += qualstr[t-i];

        for(t=l-tail-1; t-w>=front; t--) {
            totalQual += qualstr[t-w+1];
            // rolling
            if(t < l-tail-1) {
                totalQual -= qualstr[t+1];
            }
            // add 33 for phred33 transforming
            if((double)totalQual / (double)w >= 33 + mOptions->qualityCut.qualityTail)
                break;
        }

        if(t < l-1)
            t = t-w+1;
        while(t>=0 && seq[t] == 'N')
            t--;
        rlen = t - front + 1;
    }

    if(rlen <= 0 || front >= l-1)
        return NULL;

    r->mSeq->erase(0, front);
    r->mSeq->resize(rlen);
    r->mQuality->erase(0, front);
    r->mQuality->resize(rlen);

    frontTrimmed = front;

    return r;
}

bool Filter::match(vector<string>& list, string target, int threshold) {
    for(int i=0; i<list.size(); i++) {
        int diff = 0;
        int len1 = list[i].length();
        int len2 = target.length();
        for(int s=0; s<len1 && s<len2; s++) {
            if(list[i][s] != target[s]) {
                diff++;
                if(diff>threshold)
                    break;
            }
        }
        if(diff <= threshold)
            return true;
    }
    return false;
}
