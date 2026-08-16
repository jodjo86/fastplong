#include <stdlib.h>
#include "filterresult.h"
#include "stats.h"
#include "htmlreporter.h"
#include <memory.h>

FilterResult::FilterResult(Options* opt, bool paired){
    mOptions = opt;
    mTrimmedAdapterRead = 0;
    mTrimmedAdapterBases = 0;
    mSamplingDroppedReads = 0;
    mSamplingDroppedBases = 0;
    mBestReadSegmentTrimmedReads = 0;
    mBestReadSegmentTrimmedBases = 0;
    mChimericReads = 0;
    mChimeraRemovedBases = 0;
    mChimeraProducedSegments = 0;
    for(int i=0; i<FILTER_RESULT_TYPES; i++) {
        mFilterReadStats[i] = 0;
    }
    mCorrectionMatrix = new long[64];
    memset(mCorrectionMatrix, 0, sizeof(long)*64);
}

FilterResult::~FilterResult() {
    delete mCorrectionMatrix;
}

void FilterResult::addFilterResult(int result, int readNum) {
    if(result < PASS_FILTER || result >= FILTER_RESULT_TYPES)
        return ;
    mFilterReadStats[result] += readNum;
}

FilterResult* FilterResult::merge(vector<FilterResult*>& list) {
    if(list.size() == 0)
        return nullptr;
    FilterResult* result = new FilterResult(list[0]->mOptions);

    long* target = result->getFilterReadStats();
    // read stats
    for(int i=0; i<list.size(); i++) {
        long* current = list[i]->getFilterReadStats();
        for(int j=0; j<FILTER_RESULT_TYPES; j++) {
            target[j] += current[j];
        }
        result->mTrimmedAdapterRead += list[i]->mTrimmedAdapterRead;
        result->mTrimmedAdapterBases += list[i]->mTrimmedAdapterBases;
        result->mSamplingDroppedReads += list[i]->mSamplingDroppedReads;
        result->mSamplingDroppedBases += list[i]->mSamplingDroppedBases;
        result->mBestReadSegmentTrimmedReads += list[i]->mBestReadSegmentTrimmedReads;
        result->mBestReadSegmentTrimmedBases += list[i]->mBestReadSegmentTrimmedBases;
        result->mChimericReads += list[i]->mChimericReads;
        result->mChimeraRemovedBases += list[i]->mChimeraRemovedBases;
        result->mChimeraProducedSegments += list[i]->mChimeraProducedSegments;
        result->mBestReadCandidates.insert(result->mBestReadCandidates.end(), list[i]->mBestReadCandidates.begin(), list[i]->mBestReadCandidates.end());

        for(int b=0; b<4; b++) {
          result->mTrimmedPolyXReads[b] += list[i]->mTrimmedPolyXReads[b];
          result->mTrimmedPolyXBases[b] += list[i]->mTrimmedPolyXBases[b];
        }

        // merge adapter stats
        map<string, long>::iterator iter;
        for(iter = list[i]->mAdapter.begin(); iter != list[i]->mAdapter.end(); iter++) {
            if(result->mAdapter.count(iter->first) > 0)
                result->mAdapter[iter->first] += iter->second;
            else
                result->mAdapter[iter->first] = iter->second;
        }
    }

    // sort adapters list by adapter length from short to long

    return result;
}

void FilterResult::addReadTrimmed(int bases) {
    mTrimmedAdapterBases  += bases;
    mTrimmedAdapterRead++;
}

void FilterResult::addAdapterTrimmed(string adapter ) {
    if(adapter.empty())
        return;
    
    if(mAdapter.count(adapter) >0 )
        mAdapter[adapter]++;
    else
        mAdapter[adapter] = 1;
}

void FilterResult::addPolyXTrimmed(int base, int length) {
    mTrimmedPolyXReads[base] += 1;
    mTrimmedPolyXBases[base] += length;
}

void FilterResult::addSamplingDropped(int bases) {
    mSamplingDroppedReads++;
    mSamplingDroppedBases += bases;
}

void FilterResult::addBestReadCandidate(unsigned long long key, double score, int length) {
    mBestReadCandidates.push_back(BestReadRecord(key, score, length));
}

void FilterResult::addBestReadSegmentTrimmed(int bases) {
    if(bases <= 0)
        return;
    mBestReadSegmentTrimmedReads++;
    mBestReadSegmentTrimmedBases += bases;
}

void FilterResult::addChimericRead(int removedBases, int producedSegments) {
    mChimericReads++;
    if(removedBases > 0)
        mChimeraRemovedBases += removedBases;
    if(producedSegments > 0)
        mChimeraProducedSegments += producedSegments;
}

long FilterResult::getTotalPolyXTrimmedReads() {
  long sum_reads = 0;
  for(int b = 0; b < 4; b++)
    sum_reads += mTrimmedPolyXReads[b];
  return sum_reads;
}

long FilterResult::getTotalPolyXTrimmedBases() {
  long sum_bases = 0;
  for(int b = 0; b < 4; b++)
    sum_bases += mTrimmedPolyXBases[b];
  return sum_bases;
}


void FilterResult::print() {
    cerr <<  "reads passed filter: " << mFilterReadStats[PASS_FILTER] << endl;
    cerr <<  "reads failed due to low quality: " << mFilterReadStats[FAIL_QUALITY] << endl;
    cerr <<  "reads failed due to too many N: " << mFilterReadStats[FAIL_N_BASE] << endl;
    if(mOptions->gcContentFilter.enabled) {
        cerr <<  "reads failed due to GC content: " << mFilterReadStats[FAIL_GC_CONTENT] << endl;
    }
    if(mOptions->lengthFilter.enabled) {
        cerr <<  "reads failed due to too short: " << mFilterReadStats[FAIL_LENGTH] << endl;
        if(mOptions->lengthFilter.maxLength > 0)
            cerr <<  "reads failed due to too long: " << mFilterReadStats[FAIL_TOO_LONG] << endl;
    }
    if(mOptions->complexityFilter.enabled) {
        cerr <<  "reads failed due to low complexity: " << mFilterReadStats[FAIL_COMPLEXITY] << endl;
    }
    if(mOptions->sampling.enabled) {
        cerr <<  "reads discarded by sampling: " << mSamplingDroppedReads << endl;
        cerr <<  "bases discarded by sampling: " << mSamplingDroppedBases << endl;
    }
    if(mOptions->bestRead.enabled) {
        cerr <<  "candidate reads for best read selection: " << mOptions->bestRead.candidateReads << endl;
        cerr <<  "candidate bases for best read selection: " << mOptions->bestRead.candidateBases << endl;
        cerr <<  "selected reads by best read selection: " << mOptions->bestRead.selectedReads << endl;
        cerr <<  "selected bases by best read selection: " << mOptions->bestRead.selectedBases << endl;
    }
    if(mOptions->bestReadSegment.enabled) {
        cerr <<  "reads trimmed by best-read-segment mode: " << mBestReadSegmentTrimmedReads << endl;
        cerr <<  "bases trimmed by best-read-segment mode: " << mBestReadSegmentTrimmedBases << endl;
    }
    if(mOptions->adapter.enabled) {
        cerr <<  "reads with adapter trimmed: " << mTrimmedAdapterRead << endl;
        cerr <<  "bases trimmed due to adapters: " << mTrimmedAdapterBases << endl;
        if(mOptions->adapter.splitChimera) {
            cerr <<  "reads with internal adapters/chimeric signals: " << mChimericReads << endl;
            cerr <<  "bases removed by chimera splitting: " << mChimeraRemovedBases << endl;
            cerr <<  "segments produced by chimera splitting: " << mChimeraProducedSegments << endl;
            if(mOptions->adapter.discardChimera)
                cerr <<  "reads failed due to chimeric adapters: " << mFilterReadStats[FAIL_CHIMERA] << endl;
        }
    }
    if(mOptions->polyXTrim.enabled) {
        cerr <<  "reads with polyX in 3' end: " << getTotalPolyXTrimmedReads() << endl;
        cerr <<  "bases trimmed in polyX tail: " << getTotalPolyXTrimmedBases() << endl;
    }
}

void FilterResult::reportJson(ofstream& ofs, string padding) {
    ofs << "{" << endl;

    ofs << padding << "\t" << "\"passed_filter_reads\": " << mFilterReadStats[PASS_FILTER] << "," << endl;
    ofs << padding << "\t" << "\"low_quality_reads\": " << mFilterReadStats[FAIL_QUALITY] << "," << endl;
    ofs << padding << "\t" << "\"too_many_N_reads\": " << mFilterReadStats[FAIL_N_BASE] << "," << endl;
    if(mOptions->gcContentFilter.enabled)
        ofs << padding << "\t" << "\"gc_content_filtered_reads\": " << mFilterReadStats[FAIL_GC_CONTENT] << "," << endl;
    if(mOptions->complexityFilter.enabled)
        ofs << padding << "\t" << "\"low_complexity_reads\": " << mFilterReadStats[FAIL_COMPLEXITY] << "," << endl;
    ofs << padding << "\t" << "\"too_short_reads\": " << mFilterReadStats[FAIL_LENGTH] << "," << endl;
    ofs << padding << "\t" << "\"too_long_reads\": " << mFilterReadStats[FAIL_TOO_LONG];
    if(mOptions->adapter.enabled && mOptions->adapter.splitChimera) {
        ofs << "," << endl;
        ofs << padding << "\t" << "\"chimeric_reads\": " << mChimericReads << "," << endl;
        ofs << padding << "\t" << "\"chimera_removed_bases\": " << mChimeraRemovedBases << "," << endl;
        ofs << padding << "\t" << "\"chimera_produced_segments\": " << mChimeraProducedSegments;
        if(mOptions->adapter.discardChimera) {
            ofs << "," << endl;
            ofs << padding << "\t" << "\"chimeric_failed_reads\": " << mFilterReadStats[FAIL_CHIMERA];
        }
    }
    if(mOptions->sampling.enabled) {
        ofs << "," << endl;
        ofs << padding << "\t" << "\"sampling_dropped_reads\": " << mSamplingDroppedReads << "," << endl;
        ofs << padding << "\t" << "\"sampling_dropped_bases\": " << mSamplingDroppedBases << endl;
    } else if(mOptions->bestRead.enabled) {
        ofs << "," << endl;
        ofs << padding << "\t" << "\"best_selection_candidate_reads\": " << mOptions->bestRead.candidateReads << "," << endl;
        ofs << padding << "\t" << "\"best_selection_candidate_bases\": " << mOptions->bestRead.candidateBases << "," << endl;
        ofs << padding << "\t" << "\"best_selection_selected_reads\": " << mOptions->bestRead.selectedReads << "," << endl;
        ofs << padding << "\t" << "\"best_selection_selected_bases\": " << mOptions->bestRead.selectedBases << endl;
    } else if(mOptions->bestReadSegment.enabled) {
        ofs << "," << endl;
        ofs << padding << "\t" << "\"best_read_segment_trimmed_reads\": " << mBestReadSegmentTrimmedReads << "," << endl;
        ofs << padding << "\t" << "\"best_read_segment_trimmed_bases\": " << mBestReadSegmentTrimmedBases << endl;
    } else {
        ofs << endl;
    }

    ofs << padding << "}," << endl;
}

void FilterResult::outputAdaptersJson(ofstream& ofs, map<string, long, classcomp>& adapterCounts) {
    map<string, long>::iterator iter;

    long total = 0;
    for(iter = adapterCounts.begin(); iter!=adapterCounts.end(); iter++) {
        total += iter->second;
    }

    if(total == 0)
        return ;

    const double reportThreshold = 0.01;
    const double dTotal = (double)total;
    bool firstItem = true;
    long reported = 0;
    for(iter = adapterCounts.begin(); iter!=adapterCounts.end(); iter++) {
        if(iter->second /dTotal < reportThreshold )
            continue;

        if(!firstItem)
            ofs << ", ";
        else
            firstItem = false;
        ofs << "\"" << iter->first << "\":" << iter->second;

        reported += iter->second;
    }

    long unreported = total - reported;

    if(unreported > 0) {
        if(!firstItem)
            ofs << ", ";
        ofs << "\"" << "others" << "\":" << unreported;
    }
}

void FilterResult::reportAdapterJson(ofstream& ofs, string padding) {
    ofs << "{" << endl;

    ofs << padding << "\t" << "\"adapter_trimmed_reads\": " << mTrimmedAdapterRead << "," << endl;
    ofs << padding << "\t" << "\"adapter_trimmed_bases\": " << mTrimmedAdapterBases << "," << endl;
    ofs << padding << "\t" << "\"read_start_adapter\": \"" << mOptions->getReadStartAdapter() << "\"," << endl;
    ofs << padding << "\t" << "\"read_end_adapter\": \"" << mOptions->getReadEndAdapter() << "\"," << endl;

    ofs << padding << "\t" << "\"read_adapter_counts\": " << "{";
        outputAdaptersJson(ofs, mAdapter);
    ofs << "}";
    if(mOptions->adapter.splitChimera) {
        ofs << "," << endl;
        ofs << padding << "\t" << "\"chimeric_reads\": " << mChimericReads << "," << endl;
        ofs << padding << "\t" << "\"chimera_removed_bases\": " << mChimeraRemovedBases << "," << endl;
        ofs << padding << "\t" << "\"chimera_produced_segments\": " << mChimeraProducedSegments;
    }
    ofs << endl;

    ofs << padding << "}," << endl;
}

void writeBaseCountsJson(ofstream& ofs, string pad, string key, long total, long (&counts)[4]) {
  ofs << pad << "\t\"total_" << key << "\": " << total << "," << endl;
  ofs << pad << "\t\"" << key << "\":{";
  for (int b=0; b<4; b++) {
    if(b > 0)
      ofs << ", ";
    ofs << "\"" << ATCG_BASES[b] << "\": " << counts[b];
  }
  ofs << "}";
}

void FilterResult::reportPolyXTrimJson(ofstream& ofs, string padding) {
    ofs << padding << "{" << endl;
    writeBaseCountsJson(ofs, padding, "polyx_trimmed_reads", getTotalPolyXTrimmedReads(), mTrimmedPolyXReads);
    ofs << "," << endl;
    writeBaseCountsJson(ofs, padding, "polyx_trimmed_bases", getTotalPolyXTrimmedBases(), mTrimmedPolyXBases);
    ofs << endl << padding << "}," << endl;
}

/*void FilterResult::reportHtml(ofstream& ofs, long totalReads) {
    const int types = 4;
    const string divName = "filtering_result";
    string labels[4] = {"good_reads", "low_quality_reads", "too_many_N_reads", "too_short_reads"};
    long counts[4] = {mFilterReadStats[PASS_FILTER], mFilterReadStats[FAIL_QUALITY], mFilterReadStats[FAIL_N_BASE], mFilterReadStats[FAIL_LENGTH]};
    
    string json_str = "var data=[";
    json_str += "{values:[" + Stats::list2string(counts, types) + "],";
    json_str += "labels:['good_reads', 'low_quality_reads', 'too_many_N_reads', 'too_short_reads'],";
    json_str += "textinfo: 'none',";
    json_str += "type:'pie'}];\n";
    string title = "Filtering statistics of sampled " + to_string(totalReads) + " reads";
    json_str += "var layout={title:'" + title + "', width:800, height:400};\n";
    json_str += "Plotly.newPlot('" + divName + "', data, layout);\n";

    ofs << "<div class='figure' id='" + divName + "'></div>\n";
    ofs << "\n<script type=\"text/javascript\">" << endl;
    ofs << json_str;
    ofs << "</script>" << endl;
} */

void FilterResult::reportHtml(ofstream& ofs, long totalReads, long totalBases) {
    double total = (double)totalReads;
    ofs << "<table class='summary_table'>\n";
    HtmlReporter::outputRow(ofs, "reads passed filters:", HtmlReporter::formatNumber(mFilterReadStats[PASS_FILTER]) + " (" + to_string(mFilterReadStats[PASS_FILTER] * 100.0 / total) + "%)");

    HtmlReporter::outputRow(ofs, "reads with low quality:", HtmlReporter::formatNumber(mFilterReadStats[FAIL_QUALITY]) + " (" + to_string(mFilterReadStats[FAIL_QUALITY] * 100.0 / total) + "%)");
    HtmlReporter::outputRow(ofs, "reads with too many N:", HtmlReporter::formatNumber(mFilterReadStats[FAIL_N_BASE]) + " (" + to_string(mFilterReadStats[FAIL_N_BASE] * 100.0 / total) + "%)");
    if(mOptions->gcContentFilter.enabled)
        HtmlReporter::outputRow(ofs, "reads with bad GC content:", HtmlReporter::formatNumber(mFilterReadStats[FAIL_GC_CONTENT]) + " (" + to_string(mFilterReadStats[FAIL_GC_CONTENT] * 100.0 / total) + "%)");
    if(mOptions->lengthFilter.enabled) {
        HtmlReporter::outputRow(ofs, "reads too short:", HtmlReporter::formatNumber(mFilterReadStats[FAIL_LENGTH]) + " (" + to_string(mFilterReadStats[FAIL_LENGTH] * 100.0 / total) + "%)");
        if(mOptions->lengthFilter.maxLength > 0)
            HtmlReporter::outputRow(ofs, "reads too long:", HtmlReporter::formatNumber(mFilterReadStats[FAIL_TOO_LONG]) + " (" + to_string(mFilterReadStats[FAIL_TOO_LONG] * 100.0 / total) + "%)");
    }
    if(mOptions->complexityFilter.enabled)
        HtmlReporter::outputRow(ofs, "reads with low complexity:", HtmlReporter::formatNumber(mFilterReadStats[FAIL_COMPLEXITY]) + " (" + to_string(mFilterReadStats[FAIL_COMPLEXITY] * 100.0 / total) + "%)");
    if(mOptions->sampling.enabled) {
        HtmlReporter::outputRow(ofs, "reads discarded by sampling:", HtmlReporter::formatNumber(mSamplingDroppedReads) + " (" + to_string(mSamplingDroppedReads * 100.0 / total) + "%)");
        HtmlReporter::outputRow(ofs, "bases discarded by sampling:", HtmlReporter::formatNumber(mSamplingDroppedBases) + " (" + HtmlReporter::getPercents(mSamplingDroppedBases, totalBases) + "%)");
    }
    if(mOptions->bestRead.enabled) {
        HtmlReporter::outputRow(ofs, "candidate reads for best read selection:", HtmlReporter::formatNumber(mOptions->bestRead.candidateReads));
        HtmlReporter::outputRow(ofs, "candidate bases for best read selection:", HtmlReporter::formatNumber(mOptions->bestRead.candidateBases));
        HtmlReporter::outputRow(ofs, "selected reads by best read selection:", HtmlReporter::formatNumber(mOptions->bestRead.selectedReads));
        HtmlReporter::outputRow(ofs, "selected bases by best read selection:", HtmlReporter::formatNumber(mOptions->bestRead.selectedBases));
    }
    if(mOptions->bestReadSegment.enabled) {
        HtmlReporter::outputRow(ofs, "reads trimmed by best-read-segment mode:", HtmlReporter::formatNumber(mBestReadSegmentTrimmedReads));
        HtmlReporter::outputRow(ofs, "bases trimmed by best-read-segment mode:", HtmlReporter::formatNumber(mBestReadSegmentTrimmedBases) + " (" + HtmlReporter::getPercents(mBestReadSegmentTrimmedBases, totalBases) + "%)");
    }
    if(mOptions->adapter.enabled && mOptions->adapter.splitChimera && mOptions->adapter.discardChimera)
        HtmlReporter::outputRow(ofs, "reads failed due to chimeric adapters:", HtmlReporter::formatNumber(mFilterReadStats[FAIL_CHIMERA]) + " (" + to_string(mFilterReadStats[FAIL_CHIMERA] * 100.0 / total) + "%)");
    ofs << "</table>\n";
}

void FilterResult::reportAdapterHtml(ofstream& ofs, long totalBases) {
    ofs << "<div class='subsection_title' onclick=showOrHide('read1_adapters')>Adapter or bad ligation of read1</div>\n";
    ofs << "<div id='read1_adapters'>\n";
    if(mOptions->adapter.splitChimera) {
        ofs << "<table class='summary_table'>\n";
        HtmlReporter::outputRow(ofs, "reads with internal adapters/chimeric signals:", HtmlReporter::formatNumber(mChimericReads));
        HtmlReporter::outputRow(ofs, "bases removed by chimera splitting:", HtmlReporter::formatNumber(mChimeraRemovedBases) + " (" + HtmlReporter::getPercents(mChimeraRemovedBases, totalBases) + "%)");
        HtmlReporter::outputRow(ofs, "segments produced by chimera splitting:", HtmlReporter::formatNumber(mChimeraProducedSegments));
        ofs << "</table>\n";
    }
    outputAdaptersHtml(ofs, mAdapter, totalBases);
    ofs << "</div>\n";
}

void FilterResult::outputAdaptersHtml(ofstream& ofs, map<string, long, classcomp>& adapterCounts, long totalBases) {

    map<string, long>::iterator iter;

    long total = 0;
    long totalAdapterBases = 0;
    for(iter = adapterCounts.begin(); iter!=adapterCounts.end(); iter++) {
        total += iter->second;
        totalAdapterBases += iter->first.length() * iter->second;
    }

    double frac = (double)totalAdapterBases / (double)totalBases;


    if(frac < 0.01) {
        ofs << "<div class='sub_section_tips'>The input has little adapter percentage (~" << to_string(frac*100.0) << "%), probably it's trimmed before.</div>\n";
    }

    if(total == 0)
        return ;

    ofs << "<table class='summary_table'>\n";
    ofs << "<tr><td class='adapter_col' style='font-size:14px;color:#ffffff;background:#556699'>" << "Sequence" << "</td><td class='col2' style='font-size:14px;color:#ffffff;background:#556699'>" << "Occurrences" << "</td></tr>\n";

    const double reportThreshold = 0.01;
    const double dTotal = (double)total;
    long reported = 0;
    for(iter = adapterCounts.begin(); iter!=adapterCounts.end(); iter++) {
        if(iter->second /dTotal < reportThreshold )
            continue;

        ofs << "<tr><td class='adapter_col'>" << iter->first << "</td><td class='col2'>" << iter->second << "</td></tr>\n";

        reported += iter->second;
    }

    long unreported = total - reported;

    if(unreported > 0) {
        string tag = "other adapter sequences";
        if(reported == 0)
            tag = "all adapter sequences";
        ofs << "<tr><td class='adapter_col'>" << tag << "</td><td class='col2'>" << unreported << "</td></tr>\n";
    }
    ofs << "</table>\n";
}
