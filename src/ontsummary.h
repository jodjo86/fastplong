#ifndef ONT_SUMMARY_H
#define ONT_SUMMARY_H

#include <fstream>
#include <map>
#include <string>
#include <vector>

using namespace std;

class OntSummary {
public:
    OntSummary();
    bool enabled() const {return mEnabled;}
    bool parsed() const {return mParsed;}
    void setFilename(string filename);
    string filename() const {return mFilename;}
    void parse();
    void reportJson(ofstream& ofs, string padding);
    void reportHtml(ofstream& ofs);

private:
    vector<string> split(const string& line, char sep);
    int columnIndex(map<string, int>& columns, vector<string> names);
    bool parseLong(const string& s, long& value);
    bool parseDouble(const string& s, double& value);
    bool parseBool(const string& s, bool& value);
    void ensureTimeBin(int bin);
    void makeLengthBins(vector<long>& starts, vector<long>& ends, vector<long>& readCounts, vector<long>& baseCounts);
    void calcN50();
    void reportHtmlBasic(ofstream& ofs);
    void reportHtmlPassFail(ofstream& ofs);
    void reportHtmlYieldOverTime(ofstream& ofs);
    void reportHtmlReadRateAndQuality(ofstream& ofs);
    void reportHtmlLengthDistribution(ofstream& ofs);
    void reportHtmlQscoreDistribution(ofstream& ofs);
    void reportHtmlChannelActivity(ofstream& ofs);
    string list2string(vector<long>& values);
    string list2string(vector<double>& values);
    string list2string(long* values, int size);
    string list2string(double* values, int size);

private:
    bool mEnabled;
    bool mParsed;
    string mFilename;
    long mTotalReads;
    long mPassReads;
    long mFailReads;
    long mTotalBases;
    long mPassBases;
    long mFailBases;
    double mTotalQscore;
    double mMaxEndTime;
    long mN50;
    map<int, long> mLengthHist;
    long mQscoreHist[128];
    long mQscoreBases[128];
    vector<long> mTimeReads;
    vector<long> mTimeBases;
    vector<double> mTimeQscoreSum;
    vector<long> mChannelReads;
    vector<long> mChannelBases;
};

#endif
