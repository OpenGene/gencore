#ifndef GENCORE_H
#define GENCORE_H

#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "htslib/sam.h"
#include "options.h"
#include "cluster.h"
#include "stats.h"
#include <map>

using namespace std;

class Gencore {
public:
    Gencore(Options *opt);
    ~Gencore();

    void consensus();

private:
	void releaseClusters(map<int, map<int, map<int, Cluster*>>>& clusters);
	void dumpClusters(map<int, map<int, map<int, Cluster*>>>& clusters);
	void addToCluster(bam1_t* b);
	void addToProperCluster(bam1_t* b);
	void addToUnProperCluster(bam1_t* b);
	void createCluster(map<int, map<int, map<int, Cluster*>>>& clusters, int tid, int left, int right);
    void outputPair(Pair* p);
    void finishConsensus(map<int, map<int, map<int, Cluster*>>>& clusters);
    void report();

private:
    string mInput;
    string mOutput;
    Options *mOptions;
    // chrid:left:right
    map<int, map<int, map<int, Cluster*>>> mProperClusters;
    map<int, map<int, map<int, Cluster*>>> mUnProperClusters;
    bam_hdr_t *mBamHeader;
    samFile* mOutSam;
    Stats* mPreStats;
    Stats* mPostStats;
};

#endif