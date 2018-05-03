#ifndef REFERENCE_H
#define REFERENCE_H

// includes
#include "fastareader.h"
#include "options.h"

using namespace std;

// Singleton reference handler

class Reference
{
public:
    ~Reference();

    const char* getData(int contig, int pos, int len);

    static Reference* instance(Options* opt);

private:
    Reference(Options* opt);

private:
    FastaReader* mRef;
    static Reference* mInstance;
    Options* mOptions;
    int mLastBamContig;
    int mLastLen;
    const char* mLastData;
};


#endif

