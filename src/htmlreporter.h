#ifndef HTML_REPORTER_H
#define HTML_REPORTER_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "options.h"
#include "stats.h"
#include <fstream>

using namespace std;

class HtmlReporter{
public:
    HtmlReporter(Options* opt);
    ~HtmlReporter();
    void setInsertHist(long* insertHist, int insertSizePeak);
    void report( Stats* preStats, Stats* postStats);

    static void outputRow(ofstream& ofs, string key, long value);
    static void outputRow(ofstream& ofs, string key, string value);
    static void outputTripleRow(ofstream& ofs, string key, string value1, string value2);
    static string formatNumber(long number);
    static string getPercents(long numerator, long denominator);
private:
    const string getCurrentSystemTime();
    void printHeader(ofstream& ofs);
    void printCSS(ofstream& ofs);
    void printJS(ofstream& ofs);
    void printFooter(ofstream& ofs);
    void reportDuplication(ofstream& ofs, long* dupHist, Stats* preStats);
    void reportCoverage(ofstream& ofs, Stats* preStats, Stats* postStats);
    void reportInsertSize(ofstream& ofs, int isizeLimit);
    void printSummary(ofstream& ofs, Stats* preStats, Stats* postStats);
    
private:
    Options* mOptions;
    double* mDupMeanGC;
    double mDupRate;
    long* mInsertHist;
    int mInsertSizePeak;
};


#endif