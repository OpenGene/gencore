#ifndef GROUP_H
#define GROUP_H

#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "pair.h"
#include "options.h"
#include "htslib/sam.h"
#include <vector>
#include <map>
#include "stats.h"

using namespace std;

class Group {
public:
    Group(Options* opt);
    ~Group();

    void dump();
    void addPair(Pair* pair);
    void addRead(bam1_t* b);

    bool matches(Pair* p);
    Pair* consensusMerge(bool crossContig);
    bam1_t* consensusMergeBam(bool isLeft, int& diff);
    int makeConsensus(vector<bam1_t* >& reads, bam1_t* out, vector<char*>& scores, bool isLeft);


    int getLeftRef(){return mPairs[0]->getLeftRef();}
    int getRightRef(){return mPairs[0]->getRightRef();}
    int getLeftPos(){return mPairs[0]->getLeftPos();}
    int getRightPos(){return mPairs[0]->getRightPos();}
    int getTLEN(){return mPairs[0]->getTLEN();}
    Pair::MapType getMapType(){return mPairs[0]->getMapType();}
    string getUMI(){return mPairs[0]->getUMI();}

    static bool test();

private:
    static int umiDiff(const string& umi1, const string& umi2);
    static bool isDuplex(const string& umi1, const string& umi2);
    
public:
    map<string, Pair*> mPairs;
    Options* mOptions;
};

#endif