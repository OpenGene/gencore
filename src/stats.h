#ifndef STATS_H
#define STATS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include "options.h"
#include <iostream>
#include <fstream>

#define MAX_SUPPORTING_READS 100

class Stats{
public:
    Stats(Options* opt);
    ~Stats();
    void addRead(int baseNum, int mismatch, bool mapped=true);
    void addMolecule(unsigned int supportingReads, bool PE);
    void addCluster(bool hasMultiMolecule);
    void print();
    long getMappedBases();
    long getMappedReads();
    long getReads();
    long getBases();
    void reportJSON(ofstream& ofs);
    double getMappingRate();
    double getDupRate();
    long getMolecules() {return mMolecule;}
    double getMismatchRate();

public:    
	static string list2string(double* list, int size);
    static string list2string(double* list, int size, long* coords);
    static string list2string(long* list, int size);

public:
    long mReadWithMismatches;
    long mCluster;
    long mMultiMoleculeCluster;
    long mMolecule;
    long mMoleculeSE;
    long mMoleculePE;
    long* mSupportingHistgram;
	Options* mOptions;
	long mBase;
	long mBaseMismatches;
	long mBaseUnmapped;
	long mRead;
	long mReadUnmapped;
    long uncountedSupportingReads;
};

#endif