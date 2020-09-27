#ifndef GENCORE_H
#define GENCORE_H

#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "htslib/sam.h"
#include "options.h"
#include "cluster.h"
#include "stats.h"
#include "bed.h"
#include "htslib/sam.h"
#include <map>
#include <set>
#include "bamutil.h"

using namespace std;

struct bamComp{
    bool operator()(const bam1_t* b1, const bam1_t* b2) const {
        if(b1->core.tid >= 0) {        // b1 is mapped
            if(b2->core.tid<0 )
                return true;
            else if(b2->core.tid >  b1->core.tid)
                return true;
            else if(b2->core.tid == b1->core.tid && b2->core.pos >  b1->core.pos)
                return true;
            else if(b2->core.tid == b1->core.tid && b2->core.pos == b1->core.pos && b2->core.mtid >  b1->core.mtid)
                return true;
            else if(b2->core.tid == b1->core.tid && b2->core.pos == b1->core.pos && b2->core.mtid == b1->core.mtid && b2->core.mpos >  b1->core.mpos)
                return true;
            else if(b2->core.tid == b1->core.tid && b2->core.pos == b1->core.pos && b2->core.mtid == b1->core.mtid && b2->core.mpos == b1->core.mpos) {
                if(b2->core.isize > b1->core.isize)
                    return true;
                else if(b2->core.isize == b1->core.isize && (long)b2->data > (long)b1->data) return true;
                else return false;
            } else
                return false;
        } else {         // b1 is unmapped
            if(b2->core.tid<0) { // both are unmapped
                return (long)b2->data > (long)b1->data;
            }
            else
                return false;
        }
    }
};

class Gencore {
public:
    Gencore(Options *opt);
    ~Gencore();

    void consensus();

private:
	void releaseClusters(map<int, map<int, map<long, Cluster*>>>& clusters);
	void dumpClusters(map<int, map<int, map<long, Cluster*>>>& clusters);
	void addToCluster(bam1_t* b);
	void addToProperCluster(bam1_t* b);
	void addToUnProperCluster(bam1_t* b);
	void createCluster(map<int, map<int, map<long, Cluster*>>>& clusters, int tid, int left, long right);
    void outputPair(Pair* p);
    bool outputBam(bam1_t* b);
    void finishConsensus(map<int, map<int, map<long, Cluster*>>>& clusters);
    void report();
    void outputBam(bam1_t* b, bool isLeft);
    void outputOutSet();
    void writeBam(bam1_t* b);

private:
    string mInput;
    string mOutput;
    Options *mOptions;
    // chrid:left:right
    map<int, map<int, map<long, Cluster*>>> mProperClusters;
    map<int, map<int, map<long, Cluster*>>> mUnProperClusters;
    bam_hdr_t *mBamHeader;
    samFile* mOutSam;
    Stats* mPreStats;
    Stats* mPostStats;
    set<bam1_t*, bamComp> mOutSet;
    bool mOutSetCleared;
    int mProcessedTid;
    int mProcessedPos;
};

#endif
