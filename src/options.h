#ifndef OPTIONS_H
#define OPTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <set>
#include "htslib/sam.h"

using namespace std;


class Options{
public:
    Options();
    bool validate();

public:
    string input;
    string output;
    string refFile;
    string bedFile;
    string umiPrefix;
    string reportTitle;
    int maxContig;
    bam_hdr_t* bamHeader;
    bool debug;
    bool hasBedFile;
    // json file
    string jsonFile;
    // html file
    string htmlFile;

    // thresholds
    int properReadsUmiDiffThreshold;
    int unproperReadsUmiDiffThreshold;
    int clusterSizeReq;
    int baseScoreReq;
    double scorePercentReq;
    // config of quality
    int highQuality;
    int moderateQuality;
    int lowQuality;
    // scores for PE consistence check
    char scoreOfNotOverlapped;
    char scoreOfHighQualityMatch;
    char scoreOfLowQualityMatch;
    char scoreOfBothHighQualityMismatch;
    char scoreOfBothModerateQualityMismatch;
    char scoreOfBothLowQualityMismatch;
    char scoreOfUnbalancedMismatchHighQuality;
    char scoreOfUnbalancedMismatchLowQuality;

    // threshold for skipping low complexity cluster
    int skipLowComplexityClusterThreshold;

    int coverageStep;
    int bedCoverageStep;
};

#endif