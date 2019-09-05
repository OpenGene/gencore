#ifndef FASTA_READER_H
#define FASTA_READER_H

// includes
#include <cctype>
#include <clocale>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <map>
#include "options.h"

using namespace std;

class FastaReader
{
public:
    FastaReader(Options* opt, string fastaFile, bool forceUpperCase = true);
    ~FastaReader();
    bool hasNext();
    void readNext();
    void readAll();

    inline string currentID()
    {
        return mCurrentID;
    }

    inline string currentDescription()
    {
        return mCurrentDescription;
    }

    inline map<string, unsigned char*>& contigs() {
        return mAllContigs;
    }

    static unsigned char base2bits(char base);
    static char bits2base(unsigned char bits);
    static unsigned char getBase(const unsigned char* refdata, int refpos);
    static string toString(const unsigned char* refdata, int pos, int len);


public:
    unsigned char* mCurrentSequence;
    string mCurrentID ;
    string mCurrentDescription;
    long mCurrentSize;
    map<string, unsigned char*> mAllContigs;
    map<string, long> mAllContigSizes;


private:
    bool readLine();
    bool endOfLine(char c);
    void setFastaSequenceIdDescription();
    unsigned char* to4bits(const string& str);

private:
    string mFastaFile;
    ifstream mFastaFileStream;
    bool mForceUpperCase;
    Options* mOptions;
};


#endif

