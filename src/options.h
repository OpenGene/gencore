#ifndef OPTIONS_H
#define OPTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
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
    string umiPrefix;
    int clusterSizeReq;
    int maxContig;
    int properReadsUmiDiffThreshold;
    int unproperReadsUmiDiffThreshold;
    bam_hdr_t* bamHeader;
    bool debug;
};

#endif