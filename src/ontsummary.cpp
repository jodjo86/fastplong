#include "ontsummary.h"
#include "htmlreporter.h"
#include "util.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <sstream>

OntSummary::OntSummary() {
    mEnabled = false;
    mParsed = false;
    mTotalReads = 0;
    mPassReads = 0;
    mFailReads = 0;
    mTotalBases = 0;
    mPassBases = 0;
    mFailBases = 0;
    mTotalQscore = 0.0;
    mMaxEndTime = 0.0;
    mN50 = 0;
    memset(mQscoreHist, 0, sizeof(long) * 128);
    memset(mQscoreBases, 0, sizeof(long) * 128);
}

void OntSummary::setFilename(string filename) {
    mFilename = filename;
    mEnabled = !filename.empty();
}

vector<string> OntSummary::split(const string& line, char sep) {
    vector<string> fields;
    size_t start = 0;
    while(true) {
        size_t pos = line.find(sep, start);
        if(pos == string::npos) {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, pos - start));
        start = pos + 1;
    }
    return fields;
}

int OntSummary::columnIndex(map<string, int>& columns, vector<string> names) {
    for(size_t i=0; i<names.size(); i++) {
        map<string, int>::iterator iter = columns.find(names[i]);
        if(iter != columns.end())
            return iter->second;
    }
    return -1;
}

bool OntSummary::parseLong(const string& s, long& value) {
    if(s.empty())
        return false;
    char* end = NULL;
    value = strtol(s.c_str(), &end, 10);
    return end != s.c_str();
}

bool OntSummary::parseDouble(const string& s, double& value) {
    if(s.empty())
        return false;
    char* end = NULL;
    value = strtod(s.c_str(), &end);
    return end != s.c_str();
}

bool OntSummary::parseBool(const string& s, bool& value) {
    if(s == "1" || s == "True" || s == "true" || s == "TRUE") {
        value = true;
        return true;
    }
    if(s == "0" || s == "False" || s == "false" || s == "FALSE") {
        value = false;
        return true;
    }
    return false;
}

void OntSummary::ensureTimeBin(int bin) {
    if(bin < 0)
        return;
    if((int)mTimeReads.size() <= bin) {
        mTimeReads.resize(bin + 1, 0);
        mTimeBases.resize(bin + 1, 0);
        mTimeQscoreSum.resize(bin + 1, 0.0);
    }
}

void OntSummary::parse() {
    if(!mEnabled || mParsed)
        return;

    ifstream ifs(mFilename);
    if(!ifs.good())
        error_exit("failed to open ONT sequencing summary file: " + mFilename);

    string header;
    getline(ifs, header);
    if(header.empty())
        error_exit("ONT sequencing summary file is empty: " + mFilename);

    vector<string> names = split(header, '\t');
    map<string, int> columns;
    for(size_t i=0; i<names.size(); i++)
        columns[names[i]] = i;

    int lengthCol = columnIndex(columns, {"sequence_length_template", "sequence_length", "read_length", "length"});
    int qscoreCol = columnIndex(columns, {"mean_qscore_template", "mean_qscore", "qscore", "quality"});
    int passCol = columnIndex(columns, {"passes_filtering", "passed_filtering", "passes_filter"});
    int startCol = columnIndex(columns, {"start_time"});
    int durationCol = columnIndex(columns, {"duration"});
    int channelCol = columnIndex(columns, {"channel"});

    if(lengthCol < 0)
        error_exit("ONT sequencing summary requires sequence_length_template or sequence_length column");
    if(qscoreCol < 0)
        error_exit("ONT sequencing summary requires mean_qscore_template or mean_qscore column");

    string line;
    while(getline(ifs, line)) {
        if(line.empty())
            continue;
        vector<string> fields = split(line, '\t');
        if((int)fields.size() <= lengthCol || (int)fields.size() <= qscoreCol)
            continue;

        long length = 0;
        double qscore = 0.0;
        if(!parseLong(fields[lengthCol], length))
            continue;
        if(!parseDouble(fields[qscoreCol], qscore))
            continue;
        if(length < 0)
            length = 0;

        bool passed = true;
        if(passCol >= 0 && (int)fields.size() > passCol)
            parseBool(fields[passCol], passed);

        mTotalReads++;
        mTotalBases += length;
        mTotalQscore += qscore;
        if(passed) {
            mPassReads++;
            mPassBases += length;
        } else {
            mFailReads++;
            mFailBases += length;
        }

        mLengthHist[(int)length]++;
        int q = (int)(qscore + 0.5);
        q = max(0, min(127, q));
        mQscoreHist[q]++;
        mQscoreBases[q] += length;

        if(startCol >= 0 && (int)fields.size() > startCol) {
            double start = 0.0;
            if(parseDouble(fields[startCol], start)) {
                int hour = (int)floor(start / 3600.0);
                ensureTimeBin(hour);
                mTimeReads[hour]++;
                mTimeBases[hour] += length;
                mTimeQscoreSum[hour] += qscore;
                double end = start;
                if(durationCol >= 0 && (int)fields.size() > durationCol) {
                    double duration = 0.0;
                    if(parseDouble(fields[durationCol], duration))
                        end += duration;
                }
                if(end > mMaxEndTime)
                    mMaxEndTime = end;
            }
        }

        if(channelCol >= 0 && (int)fields.size() > channelCol) {
            long channel = 0;
            if(parseLong(fields[channelCol], channel) && channel >= 0) {
                if((long)mChannelReads.size() <= channel) {
                    mChannelReads.resize(channel + 1, 0);
                    mChannelBases.resize(channel + 1, 0);
                }
                mChannelReads[channel]++;
                mChannelBases[channel] += length;
            }
        }
    }

    calcN50();
    mParsed = true;
}

void OntSummary::calcN50() {
    long half = mTotalBases / 2;
    long bases = 0;
    mN50 = 0;
    for(map<int, long>::reverse_iterator iter = mLengthHist.rbegin(); iter != mLengthHist.rend(); iter++) {
        bases += (long)iter->first * iter->second;
        if(bases >= half) {
            mN50 = iter->first;
            break;
        }
    }
}

void OntSummary::makeLengthBins(vector<long>& starts, vector<long>& ends, vector<long>& readCounts, vector<long>& baseCounts) {
    if(mLengthHist.empty())
        return;
    const int maxBins = 80;
    long minLen = max(1, mLengthHist.begin()->first);
    long maxLen = max(1, mLengthHist.rbegin()->first);
    if(maxLen - minLen + 1 <= maxBins) {
        for(long len = minLen; len <= maxLen; len++) {
            starts.push_back(len);
            ends.push_back(len);
        }
    } else {
        double logMin = log((double)minLen);
        double logMax = log((double)maxLen + 1.0);
        long start = minLen;
        for(int i=0; i<maxBins && start <= maxLen; i++) {
            long end = (long)floor(exp(logMin + (logMax - logMin) * (double)(i + 1) / (double)maxBins)) - 1;
            if(i == maxBins - 1)
                end = maxLen;
            if(end < start)
                end = start;
            if(end > maxLen)
                end = maxLen;
            starts.push_back(start);
            ends.push_back(end);
            start = end + 1;
        }
    }

    readCounts.assign(starts.size(), 0);
    baseCounts.assign(starts.size(), 0);
    size_t bin = 0;
    for(map<int, long>::iterator iter = mLengthHist.begin(); iter != mLengthHist.end(); iter++) {
        long len = iter->first;
        long count = iter->second;
        while(bin < ends.size() && len > ends[bin])
            bin++;
        if(bin >= ends.size())
            break;
        if(len >= starts[bin]) {
            readCounts[bin] += count;
            baseCounts[bin] += count * len;
        }
    }
}

string OntSummary::list2string(vector<long>& values) {
    stringstream ss;
    for(size_t i=0; i<values.size(); i++) {
        ss << values[i];
        if(i + 1 < values.size())
            ss << ",";
    }
    return ss.str();
}

string OntSummary::list2string(vector<double>& values) {
    stringstream ss;
    for(size_t i=0; i<values.size(); i++) {
        ss << values[i];
        if(i + 1 < values.size())
            ss << ",";
    }
    return ss.str();
}

string OntSummary::list2string(long* values, int size) {
    stringstream ss;
    for(int i=0; i<size; i++) {
        ss << values[i];
        if(i + 1 < size)
            ss << ",";
    }
    return ss.str();
}

string OntSummary::list2string(double* values, int size) {
    stringstream ss;
    for(int i=0; i<size; i++) {
        ss << values[i];
        if(i + 1 < size)
            ss << ",";
    }
    return ss.str();
}

void OntSummary::reportJson(ofstream& ofs, string padding) {
    parse();

    vector<long> hours(mTimeReads.size(), 0);
    vector<long> cumulativeBases(mTimeBases.size(), 0);
    vector<double> meanQscoreByHour(mTimeReads.size(), 0.0);
    long bases = 0;
    for(size_t i=0; i<mTimeReads.size(); i++) {
        hours[i] = i;
        bases += mTimeBases[i];
        cumulativeBases[i] = bases;
        if(mTimeReads[i] > 0)
            meanQscoreByHour[i] = mTimeQscoreSum[i] / (double)mTimeReads[i];
    }

    vector<long> channels;
    vector<long> channelReads;
    vector<long> channelBases;
    for(size_t i=0; i<mChannelReads.size(); i++) {
        if(mChannelReads[i] == 0 && mChannelBases[i] == 0)
            continue;
        channels.push_back(i);
        channelReads.push_back(mChannelReads[i]);
        channelBases.push_back(mChannelBases[i]);
    }

    vector<long> lenStarts;
    vector<long> lenEnds;
    vector<long> lenReads;
    vector<long> lenBases;
    makeLengthBins(lenStarts, lenEnds, lenReads, lenBases);
    vector<long> lenCenters(lenStarts.size(), 0);
    for(size_t i=0; i<lenStarts.size(); i++)
        lenCenters[i] = (lenStarts[i] + lenEnds[i]) / 2;

    int qMin = 0;
    int qMax = 0;
    while(qMin < 127 && mQscoreHist[qMin] == 0)
        qMin++;
    for(int i=127; i>=0; i--) {
        if(mQscoreHist[i] > 0) {
            qMax = i;
            break;
        }
    }
    vector<long> qscore;
    vector<long> qscoreReads;
    vector<long> qscoreBases;
    if(mTotalReads > 0) {
        for(int q=qMin; q<=qMax; q++) {
            qscore.push_back(q);
            qscoreReads.push_back(mQscoreHist[q]);
            qscoreBases.push_back(mQscoreBases[q]);
        }
    }

    ofs << "{" << endl;
    ofs << padding << "\t\"file\":\"" << mFilename << "\"," << endl;
    ofs << padding << "\t\"total_reads\":" << mTotalReads << "," << endl;
    ofs << padding << "\t\"total_bases\":" << mTotalBases << "," << endl;
    ofs << padding << "\t\"pass_reads\":" << mPassReads << "," << endl;
    ofs << padding << "\t\"pass_bases\":" << mPassBases << "," << endl;
    ofs << padding << "\t\"fail_reads\":" << mFailReads << "," << endl;
    ofs << padding << "\t\"fail_bases\":" << mFailBases << "," << endl;
    ofs << padding << "\t\"mean_qscore\":" << (mTotalReads == 0 ? 0.0 : mTotalQscore / (double)mTotalReads) << "," << endl;
    ofs << padding << "\t\"n50\":" << mN50 << "," << endl;
    ofs << padding << "\t\"run_time_seconds\":" << mMaxEndTime << "," << endl;
    ofs << padding << "\t\"yield_over_time\": {" << endl;
    ofs << padding << "\t\t\"hour\":[" << list2string(hours) << "]," << endl;
    ofs << padding << "\t\t\"read_count\":[" << list2string(mTimeReads) << "]," << endl;
    ofs << padding << "\t\t\"base_count\":[" << list2string(mTimeBases) << "]," << endl;
    ofs << padding << "\t\t\"cumulative_bases\":[" << list2string(cumulativeBases) << "]," << endl;
    ofs << padding << "\t\t\"mean_qscore\":[" << list2string(meanQscoreByHour) << "]" << endl;
    ofs << padding << "\t}," << endl;
    ofs << padding << "\t\"channel_activity\": {" << endl;
    ofs << padding << "\t\t\"channel\":[" << list2string(channels) << "]," << endl;
    ofs << padding << "\t\t\"read_count\":[" << list2string(channelReads) << "]," << endl;
    ofs << padding << "\t\t\"base_count\":[" << list2string(channelBases) << "]" << endl;
    ofs << padding << "\t}," << endl;
    ofs << padding << "\t\"read_length_distribution\": {" << endl;
    ofs << padding << "\t\t\"start\":[" << list2string(lenStarts) << "]," << endl;
    ofs << padding << "\t\t\"end\":[" << list2string(lenEnds) << "]," << endl;
    ofs << padding << "\t\t\"center\":[" << list2string(lenCenters) << "]," << endl;
    ofs << padding << "\t\t\"read_count\":[" << list2string(lenReads) << "]," << endl;
    ofs << padding << "\t\t\"base_count\":[" << list2string(lenBases) << "]" << endl;
    ofs << padding << "\t}," << endl;
    ofs << padding << "\t\"qscore_distribution\": {" << endl;
    ofs << padding << "\t\t\"qscore\":[" << list2string(qscore) << "]," << endl;
    ofs << padding << "\t\t\"read_count\":[" << list2string(qscoreReads) << "]," << endl;
    ofs << padding << "\t\t\"base_count\":[" << list2string(qscoreBases) << "]" << endl;
    ofs << padding << "\t}" << endl;
    ofs << padding << "}";
}

void OntSummary::reportHtml(ofstream& ofs) {
    parse();
    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('ont_summary')><a name='summary'>ONT sequencing summary</a></div>\n";
    ofs << "<div id='ont_summary'>\n";
    reportHtmlBasic(ofs);
    reportHtmlPassFail(ofs);
    reportHtmlYieldOverTime(ofs);
    reportHtmlReadRateAndQuality(ofs);
    reportHtmlLengthDistribution(ofs);
    reportHtmlQscoreDistribution(ofs);
    reportHtmlChannelActivity(ofs);
    ofs << "</div>\n";
    ofs << "</div>\n";
}

void OntSummary::reportHtmlBasic(ofstream& ofs) {
    ofs << "<div class='subsection_title'>ONT summary: Basic statistics</div>\n";
    ofs << "<table>\n";
    HtmlReporter::outputRow(ofs, "sequencing summary:", mFilename);
    HtmlReporter::outputRow(ofs, "total reads:", HtmlReporter::formatNumber(mTotalReads));
    HtmlReporter::outputRow(ofs, "total bases:", HtmlReporter::formatNumber(mTotalBases));
    HtmlReporter::outputRow(ofs, "pass reads:", HtmlReporter::formatNumber(mPassReads) + " (" + HtmlReporter::getPercents(mPassReads, mTotalReads) + "%)");
    HtmlReporter::outputRow(ofs, "pass bases:", HtmlReporter::formatNumber(mPassBases) + " (" + HtmlReporter::getPercents(mPassBases, mTotalBases) + "%)");
    HtmlReporter::outputRow(ofs, "mean qscore:", to_string(mTotalReads == 0 ? 0.0 : mTotalQscore / (double)mTotalReads));
    HtmlReporter::outputRow(ofs, "N50 length:", HtmlReporter::formatNumber(mN50));
    HtmlReporter::outputRow(ofs, "run time:", to_string(mMaxEndTime / 3600.0) + " hours");
    ofs << "</table>\n";
}

void OntSummary::reportHtmlPassFail(ofstream& ofs) {
    ofs << "<div class='subsection_title'>ONT summary: Pass/fail yield</div>\n";
    ofs << "<div class='figure' id='plot_ont_pass_fail' style='height:420px;'></div>\n";
    ofs << "\n<script type=\"text/javascript\">" << endl;
    ofs << "var ontPassFailReads={x:['pass','fail'],y:[" << mPassReads << "," << mFailReads << "],name:'reads',type:'bar',marker:{color:'rgb(73,120,180)'}};\n";
    ofs << "var ontPassFailBases={x:['pass','fail'],y:[" << mPassBases << "," << mFailBases << "],name:'bases',type:'bar',yaxis:'y2',marker:{color:'rgb(77,156,96)'}};\n";
    ofs << "var layout={title:'ONT pass/fail yield',barmode:'group',yaxis:{title:'reads'},yaxis2:{title:'bases',overlaying:'y',side:'right'}};\n";
    ofs << "Plotly.newPlot('plot_ont_pass_fail',[ontPassFailReads,ontPassFailBases],layout);\n";
    ofs << "</script>" << endl;
}

void OntSummary::reportHtmlYieldOverTime(ofstream& ofs) {
    vector<long> hours(mTimeReads.size(), 0);
    vector<long> cumulativeBases(mTimeBases.size(), 0);
    long bases = 0;
    for(size_t i=0; i<mTimeReads.size(); i++) {
        hours[i] = i;
        bases += mTimeBases[i];
        cumulativeBases[i] = bases;
    }
    ofs << "<div class='subsection_title'>ONT summary: Yield over time</div>\n";
    ofs << "<div class='figure' id='plot_ont_yield_time' style='height:420px;'></div>\n";
    ofs << "\n<script type=\"text/javascript\">" << endl;
    ofs << "var ontHourlyBases={x:[" << list2string(hours) << "],y:[" << list2string(mTimeBases) << "],name:'bases per hour',type:'bar',marker:{color:'rgb(77,156,96)'}};\n";
    ofs << "var ontCumulativeBases={x:[" << list2string(hours) << "],y:[" << list2string(cumulativeBases) << "],name:'cumulative bases',type:'scatter',mode:'lines',yaxis:'y2',line:{color:'rgb(192,105,77)',width:2}};\n";
    ofs << "var layout={title:'ONT yield over time',xaxis:{title:'run time (hour)'},yaxis:{title:'bases per hour'},yaxis2:{title:'cumulative bases',overlaying:'y',side:'right'}};\n";
    ofs << "Plotly.newPlot('plot_ont_yield_time',[ontHourlyBases,ontCumulativeBases],layout);\n";
    ofs << "</script>" << endl;
}

void OntSummary::reportHtmlReadRateAndQuality(ofstream& ofs) {
    vector<long> hours(mTimeReads.size(), 0);
    vector<double> meanQscore(mTimeReads.size(), 0.0);
    for(size_t i=0; i<mTimeReads.size(); i++) {
        hours[i] = i;
        if(mTimeReads[i] > 0)
            meanQscore[i] = mTimeQscoreSum[i] / (double)mTimeReads[i];
    }
    ofs << "<div class='subsection_title'>ONT summary: Read count and quality over time</div>\n";
    ofs << "<div class='figure' id='plot_ont_reads_quality_time' style='height:420px;'></div>\n";
    ofs << "\n<script type=\"text/javascript\">" << endl;
    ofs << "var ontHourlyReads={x:[" << list2string(hours) << "],y:[" << list2string(mTimeReads) << "],name:'reads per hour',type:'bar',marker:{color:'rgb(73,120,180)'}};\n";
    ofs << "var ontHourlyQscore={x:[" << list2string(hours) << "],y:[" << list2string(meanQscore) << "],name:'mean qscore',type:'scatter',mode:'lines+markers',yaxis:'y2',line:{color:'rgb(102,82,163)',width:2},marker:{color:'rgb(102,82,163)',size:5}};\n";
    ofs << "var layout={title:'ONT read count and quality over time',xaxis:{title:'run time (hour)'},yaxis:{title:'reads per hour'},yaxis2:{title:'mean qscore',overlaying:'y',side:'right'}};\n";
    ofs << "Plotly.newPlot('plot_ont_reads_quality_time',[ontHourlyReads,ontHourlyQscore],layout);\n";
    ofs << "</script>" << endl;
}

void OntSummary::reportHtmlLengthDistribution(ofstream& ofs) {
    vector<long> starts;
    vector<long> ends;
    vector<long> reads;
    vector<long> bases;
    makeLengthBins(starts, ends, reads, bases);
    vector<long> centers(starts.size(), 0);
    vector<double> readPercents(starts.size(), 0.0);
    vector<double> basePercents(starts.size(), 0.0);
    for(size_t i=0; i<starts.size(); i++) {
        centers[i] = (starts[i] + ends[i]) / 2;
        readPercents[i] = mTotalReads == 0 ? 0.0 : (double)reads[i] * 100.0 / (double)mTotalReads;
        basePercents[i] = mTotalBases == 0 ? 0.0 : (double)bases[i] * 100.0 / (double)mTotalBases;
    }
    ofs << "<div class='subsection_title'>ONT summary: Read length distribution</div>\n";
    ofs << "<div class='figure' id='plot_ont_length_distribution' style='height:420px;'></div>\n";
    ofs << "\n<script type=\"text/javascript\">" << endl;
    ofs << "var ontLengthReads={x:[" << list2string(centers) << "],y:[" << list2string(readPercents) << "],name:'% reads',type:'scatter',mode:'lines+markers',line:{color:'rgb(73,120,180)',width:2},marker:{color:'rgb(73,120,180)',size:5}};\n";
    ofs << "var ontLengthBases={x:[" << list2string(centers) << "],y:[" << list2string(basePercents) << "],name:'% bases',type:'scatter',mode:'lines+markers',line:{color:'rgb(77,156,96)',width:2},marker:{color:'rgb(77,156,96)',size:5}};\n";
    ofs << "var layout={title:'ONT read length distribution',xaxis:{title:'read length',type:'log'},yaxis:{title:'Percent (%)',rangemode:'tozero'}};\n";
    ofs << "Plotly.newPlot('plot_ont_length_distribution',[ontLengthReads,ontLengthBases],layout);\n";
    ofs << "</script>" << endl;
}

void OntSummary::reportHtmlQscoreDistribution(ofstream& ofs) {
    int qMin = 0;
    int qMax = 0;
    while(qMin < 127 && mQscoreHist[qMin] == 0)
        qMin++;
    for(int i=127; i>=0; i--) {
        if(mQscoreHist[i] > 0) {
            qMax = i;
            break;
        }
    }
    vector<long> qscores;
    vector<double> readPercents;
    vector<double> basePercents;
    if(mTotalReads > 0) {
        for(int q=qMin; q<=qMax; q++) {
            qscores.push_back(q);
            readPercents.push_back((double)mQscoreHist[q] * 100.0 / (double)mTotalReads);
            basePercents.push_back(mTotalBases == 0 ? 0.0 : (double)mQscoreBases[q] * 100.0 / (double)mTotalBases);
        }
    }
    ofs << "<div class='subsection_title'>ONT summary: Read mean qscore distribution</div>\n";
    ofs << "<div class='figure' id='plot_ont_qscore_distribution' style='height:420px;'></div>\n";
    ofs << "\n<script type=\"text/javascript\">" << endl;
    ofs << "var ontQscoreReads={x:[" << list2string(qscores) << "],y:[" << list2string(readPercents) << "],name:'% reads',type:'bar',marker:{color:'rgb(73,120,180)'}};\n";
    ofs << "var ontQscoreBases={x:[" << list2string(qscores) << "],y:[" << list2string(basePercents) << "],name:'% bases',type:'bar',marker:{color:'rgb(102,82,163)'}};\n";
    ofs << "var layout={title:'ONT read mean qscore distribution',barmode:'group',xaxis:{title:'mean qscore'},yaxis:{title:'Percent (%)'}};\n";
    ofs << "Plotly.newPlot('plot_ont_qscore_distribution',[ontQscoreReads,ontQscoreBases],layout);\n";
    ofs << "</script>" << endl;
}

void OntSummary::reportHtmlChannelActivity(ofstream& ofs) {
    vector<long> channels;
    vector<long> reads;
    vector<long> bases;
    for(size_t i=0; i<mChannelReads.size(); i++) {
        if(mChannelReads[i] == 0 && mChannelBases[i] == 0)
            continue;
        channels.push_back(i);
        reads.push_back(mChannelReads[i]);
        bases.push_back(mChannelBases[i]);
    }
    ofs << "<div class='subsection_title'>ONT summary: Channel activity</div>\n";
    ofs << "<div class='figure' id='plot_ont_channel_activity' style='height:420px;'></div>\n";
    ofs << "\n<script type=\"text/javascript\">" << endl;
    ofs << "var ontChannelReads={x:[" << list2string(channels) << "],y:[" << list2string(reads) << "],name:'reads',type:'scatter',mode:'markers',marker:{color:'rgb(73,120,180)',size:5}};\n";
    ofs << "var ontChannelBases={x:[" << list2string(channels) << "],y:[" << list2string(bases) << "],name:'bases',type:'scatter',mode:'markers',yaxis:'y2',marker:{color:'rgb(77,156,96)',size:5}};\n";
    ofs << "var layout={title:'ONT channel activity',xaxis:{title:'channel'},yaxis:{title:'reads'},yaxis2:{title:'bases',overlaying:'y',side:'right'}};\n";
    ofs << "Plotly.newPlot('plot_ont_channel_activity',[ontChannelReads,ontChannelBases],layout);\n";
    ofs << "</script>" << endl;
}
