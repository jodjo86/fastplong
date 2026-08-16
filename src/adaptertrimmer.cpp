#include "adaptertrimmer.h"
#include "editdistance.h"
#include <math.h>
#include "hwy/highway.h"

AdapterTrimmer::AdapterTrimmer(){
}


AdapterTrimmer::~AdapterTrimmer(){
}

bool AdapterTrimmer::findMiddleAdapters(Read* r, string& startAdater, string& endAdapter, int& start, int& len, double edMax, int trimmingExtension) {
    len = -1;

    int startAdaterPos = searchAdapter(r->mSeq, startAdater, edMax);
    int endAdapterPos = searchAdapter(r->mSeq, endAdapter, edMax);

    if(startAdaterPos >=0 && endAdapterPos>=0) {
        start = min(startAdaterPos, endAdapterPos);
        int end = max(startAdaterPos + startAdater.length(), endAdapterPos + endAdapter.length());

        start = max(0, start-trimmingExtension);
        end = min(r->length(), end+trimmingExtension);
        len = end - start;
        return true;
    } if(startAdaterPos >=0){
        int end = min(r->length(), startAdaterPos + (int)startAdater.length() +trimmingExtension);
        start = max(0, startAdaterPos-trimmingExtension);
        len = end - start;
        return true;
    } if(endAdapterPos >=0){
        int end = min(r->length(), endAdapterPos + (int)endAdapter.length() +trimmingExtension);
        start = max(0, endAdapterPos-trimmingExtension);
        len = end - start;
        return true;
    }

    return false;
}

bool AdapterTrimmer::findNextMiddleAdapter(Read* r, vector<string>& adapters, int searchStart, int& adapterStart, int& adapterLength, double edMax) {
    adapterStart = -1;
    adapterLength = 0;
    int readLen = r->length();
    for(size_t i=0; i<adapters.size(); i++) {
        if(adapters[i].empty())
            continue;
        int alen = adapters[i].length();
        if(searchStart + alen > readLen)
            continue;
        int pos = searchAdapter(r->mSeq, adapters[i], edMax, searchStart, readLen - searchStart, true, false);
        if(pos < 0)
            continue;
        if(adapterStart < 0 || pos < adapterStart || (pos == adapterStart && alen > adapterLength)) {
            adapterStart = pos;
            adapterLength = alen;
        }
    }
    return adapterStart >= 0;
}

vector<Read*> AdapterTrimmer::splitByMiddleAdapters(Read* r, FilterResult* fr, vector<string>& adapters, double edMax, int trimmingExtension, int minSegmentLength, int* removedBases) {
    vector<Read*> out;
    if(removedBases)
        *removedBases = 0;
    if(r == NULL || adapters.empty())
        return out;

    int readLen = r->length();
    vector<pair<int, int>> adapterRegions;
    int searchStart = 0;
    while(searchStart < readLen) {
        int adapterStart = -1;
        int adapterLength = 0;
        if(!findNextMiddleAdapter(r, adapters, searchStart, adapterStart, adapterLength, edMax))
            break;
        int regionStart = max(0, adapterStart - trimmingExtension);
        int regionEnd = min(readLen - 1, adapterStart + adapterLength + trimmingExtension - 1);
        if(adapterRegions.empty() || regionStart > adapterRegions.back().second + 1) {
            adapterRegions.push_back(make_pair(regionStart, regionEnd));
        } else {
            adapterRegions.back().second = max(adapterRegions.back().second, regionEnd);
        }
        searchStart = adapterStart + max(1, adapterLength);
    }

    if(adapterRegions.empty()) {
        out.push_back(r);
        return out;
    }

    int removed = 0;
    int segmentIndex = 1;
    int lastEnd = -1;
    for(size_t i=0; i<adapterRegions.size(); i++) {
        int start = max(0, adapterRegions[i].first);
        int end = min(readLen - 1, adapterRegions[i].second);
        if(start > lastEnd + 1) {
            int len = start - lastEnd - 1;
            if(len >= minSegmentLength) {
                string* seq = new string(*r->mSeq, lastEnd + 1, len);
                string* qual = new string(*r->mQuality, lastEnd + 1, len);
                string* name = new string(*r->mName);
                name->insert(1, "chimera-segment-" + to_string(segmentIndex) + "-");
                string* strand = new string(*r->mStrand);
                out.push_back(new Read(name, seq, strand, qual));
                segmentIndex++;
            } else {
                removed += len;
            }
        }
        removed += end - start + 1;
        lastEnd = max(lastEnd, end);
    }
    if(lastEnd < readLen - 1) {
        int len = readLen - lastEnd - 1;
        if(len >= minSegmentLength) {
            string* seq = new string(*r->mSeq, lastEnd + 1, len);
            string* qual = new string(*r->mQuality, lastEnd + 1, len);
            string* name = new string(*r->mName);
            name->insert(1, "chimera-segment-" + to_string(segmentIndex) + "-");
            string* strand = new string(*r->mStrand);
            out.push_back(new Read(name, seq, strand, qual));
        } else {
            removed += len;
        }
    }

    if(removedBases)
        *removedBases = removed;
    if(fr)
        fr->addChimericRead(removed, out.size());
    return out;
}

int AdapterTrimmer::trimByMultiSequences(Read* r, FilterResult* fr, vector<string>& adapterList, double edMax, int trimmingExtension) {
    int matchReq = 4;
    if(adapterList.size() > 16)
        matchReq = 5;
    if(adapterList.size() > 256)
        matchReq = 6;
    int trimmed = 0;

    string* originalSeq = r->mSeq;
    for(int i=0; i<adapterList.size(); i++) {
        trimmed += trimBySequenceStart(r, fr, adapterList[i], edMax, trimmingExtension);
        trimmed += trimBySequenceEnd(r, fr, adapterList[i], edMax, trimmingExtension);
    }

    return trimmed;
}

int AdapterTrimmer::searchAdapter(string *read, const string &adapter, double edMax, int searchStart, int searchLen, bool asLeftAsPossible, bool asRightAsPossible)
{
    if(adapter.empty())
        return -1;

    namespace hn = hwy::HWY_NAMESPACE;
    hn::ScalableTag<uint8_t> d8;
    const size_t N = hn::Lanes(d8);

    int minMismatch = 99999; // initialized with a large mismatch
    int pos = -1;
    // for the best match
    const char *adata = adapter.c_str();
    const char *rdata = read->c_str();
    int rlen = read->length();
    int alen = adapter.length();

    int threshold = round(edMax * alen);

    int searchEnd = rlen;
    if (searchLen > 0)
    {
        searchEnd = min(rlen, searchLen + searchStart);
    }

    if (searchStart + alen > rlen)
        return -1;

    if (asLeftAsPossible)
    {
        // go from left, and return immediatedly if find a mismatch < threshold
        for (int p = searchStart; p < searchEnd - alen; p++)
        {
            size_t mismatch = 0;
            for (size_t i = 0; i < alen; i += N)
            {
                const size_t lanesToLoad = min(alen - i, N);
                const auto rdata_v = hn::LoadN(d8, reinterpret_cast<const uint8_t *>(rdata + i + p), lanesToLoad);
                const auto adata_v = hn::LoadN(d8, reinterpret_cast<const uint8_t *>(adata + i), lanesToLoad);
                const auto mismatch_mask = rdata_v != adata_v;
                mismatch += hn::CountTrue(d8, mismatch_mask);
            }

            if (mismatch <= threshold)
            {
                return p;
            }
            if (mismatch <= minMismatch)
            {
                minMismatch = mismatch;
                pos = p;
            }
        }
    }
    else if (asRightAsPossible && searchEnd > alen)
    {
        // go from right, and return immediatedly if find a mismatch <= threshold
        for (int p = searchEnd - alen; p >= searchStart; p--)
        {
            size_t mismatch = 0;
            for (size_t i = 0; i < alen; i += N)
            {
                const size_t lanesToLoad = min(alen - i, N);
                const auto rdata_v = hn::LoadN(d8, reinterpret_cast<const uint8_t *>(rdata + i + p), lanesToLoad);
                const auto adata_v = hn::LoadN(d8, reinterpret_cast<const uint8_t *>(adata + i), lanesToLoad);
                const auto mismatch_mask = rdata_v != adata_v;
                mismatch += hn::CountTrue(d8, mismatch_mask);
            }
            if (mismatch <= threshold)
            {
                return p;
            }
            if (mismatch <= minMismatch)
            {
                minMismatch = mismatch;
                pos = p;
            }
        }
    }
    else
    {
        for (int p = searchStart; p < searchEnd - alen; p++)
        {
            size_t mismatch = 0;
            for (size_t i = 0; i < alen; i += N)
            {
                const size_t lanesToLoad = min(alen - i, N);
                const auto rdata_v = hn::LoadN(d8, reinterpret_cast<const uint8_t *>(rdata + i + p), lanesToLoad);
                const auto adata_v = hn::LoadN(d8, reinterpret_cast<const uint8_t *>(adata + i), lanesToLoad);
                const auto mismatch_mask = rdata_v != adata_v;
                mismatch += hn::CountTrue(d8, mismatch_mask);
            }
            if (mismatch < minMismatch)
            {
                minMismatch = mismatch;
                pos = p;
            }
        }
    }

    if (pos >= 0)
    {
        int ed = edit_distance(rdata + pos, alen, adata, alen);
        if (ed <= threshold)
            return pos;
        else
            return -1;
    }
    else
        return -1;
}

int AdapterTrimmer::trimBySequenceStart(Read* r, FilterResult* fr, string& adapterseq, double edMax, int trimmingExtension) {
    const int WINDOW = 200;
    const int PATTERN_LEN = 16;

    int rlen = r->length();
    int alen = adapterseq.length();

    const char* adata = adapterseq.c_str();
    const char* rdata = r->mSeq->c_str();

    if(rlen < PATTERN_LEN)
        return 0;

    int plen = min(PATTERN_LEN, alen);

    // search by full match
    int mpos = searchAdapter(r->mSeq, adapterseq, edMax, 0, WINDOW, false, true);
    if(mpos >= 0) {
        // extend to make a cleaner trimming
        mpos = min(mpos + trimmingExtension, rlen - alen);
        if(fr)
            fr->addAdapterTrimmed(adapterseq);
        r->trimFront(mpos + alen);
        //cout << "L " << pos << endl;
        return mpos+alen;
    }

    // adapter not found by above full match
    // search part adapter

    int mined = -1;
    // reset pos;
    int pos = -1;
    //from tail to front, search by partly match of edit distance
    for(int p=0; p<rlen-plen && p<WINDOW - plen; p++) {
        int ed = edit_distance(rdata + p, plen, adata + alen - plen, plen);
        if(ed <= round(edMax * plen)) {
            if(pos < 0) {
                pos = p;
                mined = ed;
            } else if(ed >= mined) { // last one is best
                //break;
            }
            else {
                pos = p;
                mined = ed;
            }
        }
    }

    if(pos >= 0) {
        // extend to compare the whole adapter
        int cmplen = min(pos+plen, alen);
        int ed = edit_distance(rdata + pos + plen - cmplen, cmplen, adata + alen - cmplen, cmplen);
        if( ed <= round(edMax * cmplen) ){
            // extend to make a cleaner trimming
             pos = min(pos + trimmingExtension, rlen - alen);
            if(fr)
                fr->addAdapterTrimmed(adapterseq.substr(alen - cmplen, cmplen));
            //cout << r->mSeq->substr(0, pos + plen) << endl;
            //cout << *r->mName << endl;
            r->trimFront(pos + plen);
            //cout << "S " << pos << ", ed: " << ed << endl;
            return pos+plen;
        }
    }

    return 0;
}

int AdapterTrimmer::trimBySequenceEnd(Read* r, FilterResult* fr, string& adapterseq, double edMax, int trimmingExtension) {
    const int WINDOW = 200;
    const int PATTERN_LEN = 16;

    int rlen = r->length();
    int alen = adapterseq.length();

    const char* adata = adapterseq.c_str();
    const char* rdata = r->mSeq->c_str();

    if(rlen < PATTERN_LEN)
        return false;

    int plen = min(PATTERN_LEN, alen);

    // search by full match
    int searchStart = max(0, rlen - WINDOW);
    int mpos = searchAdapter(r->mSeq, adapterseq, edMax, searchStart, WINDOW, true, false);
    if(mpos >= 0) {
        // extend to make a cleaner trimming
        mpos = max(0, mpos - trimmingExtension);
        if(fr)
            fr->addAdapterTrimmed(adapterseq);
        r->resize(mpos);
        //cout << "R " << pos << endl;
        return rlen - mpos;
    }

    // adapter not found by above full match
    // search part adapter

    int mined = -1;
    // reset pos;
    int pos = -1;
    //from tail to front
    for(int p=0; p<rlen-plen && p<WINDOW - plen; p++) {
        int ed = edit_distance(rdata + rlen - plen -p, plen, adata, plen);
        if(ed <= round(edMax * plen)) {
            if(pos < 0) {
                pos = p;
                mined = ed;
            } else if(ed > mined) // last one is best
                break;
            else {
                pos = p;
                mined = ed;
            }
        }
    }

    if(pos > 0) {
        // extend to compare the whole adapter
        int cmplen = min(pos+plen, alen);
        if(edit_distance(rdata + rlen -plen - pos, cmplen, adata, cmplen) <= round(edMax * cmplen) ){
            // extend to make a cleaner trimming
            pos = min(pos + trimmingExtension, rlen - plen);
            if(fr)
                fr->addAdapterTrimmed(adapterseq.substr(0, cmplen));
            r->resize(rlen - plen -pos);
            //cout << "E " << pos << endl;
            return pos + plen;
        }
    }
    return 0;
}
