#ifndef PAIR_H
#define PAIR_H

#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "htslib/sam.h"
#include "options.h"

using namespace std;

class Pair {
public:
    Pair(Options* opt);
    ~Pair();

    enum MapType{Unknown, ProperlyMapped, CrossRefMapped, OnlyLeftMapped, OnlyRightMapped, NoneMapped};

    int getLeftRef();
    int getRightRef();
    int getLeftPos();
    int getRightPos();
    int getTLEN();

    char* getLeftScore();
    char* getRightScore();

    void setLeft(bam1_t *b);
    void setRight(bam1_t *b);
    bool pairFound();
    MapType getMapType();
    string getUMI();
    string getQName();
    string getLeftCigar();
    string getRightCigar();

    bool isDupWith(Pair* other);
    
    void dump();

private:
    void computeScore();

public:
    bam1_t *mLeft;
    bam1_t *mRight;
    int mMergeReads;
    int mMergeLeftDiff;
    int mMergeRightDiff;

private:
    int mTLEN;
    MapType mMapType;
    string mUMI;
    string mLeftCigar;
    string mRightCigar;
    Options* mOptions;
    char* mLeftScore;
    char* mRightScore;
};

#endif