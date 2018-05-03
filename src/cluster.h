#ifndef CLUSTER_H
#define CLUSTER_H

#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "pair.h"
#include "options.h"
#include "htslib/sam.h"
#include <vector>
#include <map>

using namespace std;

class Cluster {
public:
    Cluster(Options* opt);
    ~Cluster();

    void dump();
    void addPair(Pair* pair);
    void addRead(bam1_t* b);

    bool matches(Pair* p);
    vector<Pair*> clusterByUMI(int umiDiffThreshold = 1);
    Pair* consensusMerge();
    bam1_t* consensusMergeBam(bool isLeft, int& diff);
    int makeConsensus(vector<bam1_t* >& reads, bam1_t* out, bool isLeft);


    int getLeftRef(){return mPairs[0]->getLeftRef();}
    int getRightRef(){return mPairs[0]->getRightRef();}
    int getLeftPos(){return mPairs[0]->getLeftPos();}
    int getRightPos(){return mPairs[0]->getRightPos();}
    int getTLEN(){return mPairs[0]->getTLEN();}
    Pair::MapType getMapType(){return mPairs[0]->getMapType();}
    string getUMI(){return mPairs[0]->getUMI();}
    
public:
    vector<Pair*> mPairs;
    Options* mOptions;
};

#endif