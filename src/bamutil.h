#ifndef BAM_UTIL_H
#define BAM_UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "htslib/sam.h"

using namespace std;

class BamUtil {
public:
    BamUtil();
    ~BamUtil();

public:
    static string getQName(bam1_t *b);
    static string getUMI(string qname, const string& prefix);
    static string getUMI(bam1_t *b, const string& prefix);
    static string getSeq(bam1_t *b);
    static string getQual(bam1_t *b);
    static string getCigar(bam1_t *b);
    static char fourbits2base(uint8_t val);
    static uint8_t base2fourbits(char base);
    static void dump(bam1_t *b);
    static bool isPartOf(bam1_t *part, bam1_t *whole, bool isLeft);
    static void dumpHeader(bam_hdr_t* hdr);
    static int getRefOffset(bam1_t *b, int bampos);
    static void copyQName(bam1_t *from, bam1_t *to);
    static bool isPrimary(bam1_t *b);
    static bool isProperPair(bam1_t *b);
    static int getRightRefPos(bam1_t *b);
    static void getMOffsetAndLen(bam1_t *b, int& MOffset, int& MLen);

    static bool test();

};

#endif