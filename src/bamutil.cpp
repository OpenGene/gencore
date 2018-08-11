#include "bamutil.h"
#include <sstream>
#include <vector>
#include <memory.h>

BamUtil::BamUtil(){
}

BamUtil::~BamUtil(){
}

void BamUtil::dump(bam1_t *b) {
    cerr << b->core.tid << ":" << b->core.pos << " TLEN:" << b->core.isize << endl;
    cerr << getQName(b) << " " << getCigar(b) << endl;
    cerr << getSeq(b) << endl;
    cerr << getQual(b) << endl;
}

string BamUtil::getQName(bam1_t *b) {
    return string(bam_get_qname(b));
}

string BamUtil::getUMI(bam1_t *b, const string& prefix) {
    return getUMI(string(bam_get_qname(b)), prefix);
}

string BamUtil::getUMI(string qname, const string& prefix) {
    int len = qname.length();
    int sep = len-1;
    int prefixLen = prefix.length();

    bool foundSep = false;
    bool found= false;
    for(sep = len-1; sep>=0; sep--) {
        char c = qname[sep];
        if(c == ':') {
            foundSep = true;
            break;
        }
    }

    if(!foundSep || sep + prefixLen >=len-1)
        return "";


    // check prefix
    bool goodPrefix = true;
    for(int p=0; p<prefixLen; p++) {
        if(prefix[p] != qname[sep + 1 + p]) {
            goodPrefix = false;
            break;
        }
    }

    if(!goodPrefix)
        return "";

    int start = sep + 1 + prefixLen;
    if(start < len-1 && qname[start] == '_')
        start++;

    int numOfUnderscore = 0;
    for(int i=start; i<len; i++) {
        char c = qname[i];
        // UMI can be only A/T/C/G/N/_
        if(c != 'A' && c != 'T' && c != 'C' && c != 'G' && c != '_') {
            return "";
        }
        if(c == '_') {
            numOfUnderscore++;
            if(numOfUnderscore > 1)
                return "";
        }
    }
    return qname.substr(start, len-start);
}

string BamUtil::getQual(bam1_t *b) {
    uint8_t *data = bam_get_qual(b);
    int len = b->core.l_qseq;
    string s(len, '\0');
    for(int i=0; i<len; i++) {
        s[i] = (char)(data[i] + 33);
    }
    return s;
}

string BamUtil::getSeq(bam1_t *b) {
    uint8_t *data = bam_get_seq(b);
    int len = b->core.l_qseq;
    string s(len, '\0');
    for(int i=0; i<len; i++) {
        char base;
        if(i%2 == 1)
            base = fourbits2base(data[i/2] & 0xF);
        else
            base = fourbits2base((data[i/2]>>4) & 0xF);
        s[i] = base;
    }
    return s;
}

//Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G, 8 for T and 15 for N.
char BamUtil::fourbits2base(uint8_t val) {
    switch(val) {
        case 1:
            return 'A';
        case 2:
            return 'C';
        case 4:
            return 'G';
        case 8:
            return 'T';
        case 15:
            return 'N';
        default:
            cerr << "ERROR: Wrong base with value "<< (int)val << endl ;
            return 'N';
    }
}

uint8_t BamUtil::base2fourbits(char base) {
    switch(base) {
        case 'A':
            return 1;
        case 'C':
            return 2;
        case 'G':
            return 4;
        case 'T':
            return 8;
        case 'N':
            return 15;
        default:
            cerr << "ERROR: Wrong base "<< base << endl ;
            return 15;
    }
}

/*
@discussion In the CIGAR array, each element is a 32-bit integer. The
 lower 4 bits gives a CIGAR operation and the higher 28 bits keep the
 length of a CIGAR.
*/

string BamUtil::getCigar(bam1_t *b) {
    uint32_t *data = (uint32_t *)bam_get_cigar(b);
    int cigarNum = b->core.n_cigar;
    stringstream ss;
    for(int i=0; i<cigarNum; i++) {
        uint32_t val = data[i];
        char op = bam_cigar_opchr(val);
        uint32_t len = bam_cigar_oplen(val);
        ss << op << len;
    }
    return ss.str();
}

bool BamUtil::isPartOf(bam1_t *part, bam1_t *whole, bool isLeft) {
    uint32_t *cigarPart = (uint32_t *)bam_get_cigar(part);
    int cigarNumPart = part->core.n_cigar;
    uint32_t *cigarWhole = (uint32_t *)bam_get_cigar(whole);
    int cigarNumWhole = whole->core.n_cigar;

    if(cigarNumWhole < cigarNumPart)
        return false;

    for(int i=0; i<cigarNumPart; i++) {
        uint32_t valPart = cigarPart[i];
        // if it's right aligned, backward
        if(!isLeft)
            valPart = cigarPart[cigarNumPart - i - 1];
        char opPart = bam_cigar_op(valPart);
        uint32_t lenPart = bam_cigar_oplen(valPart);

        uint32_t valWhole = cigarWhole[i];
        // if it's right aligned, backward
        if(!isLeft)
            valWhole = cigarWhole[cigarNumWhole - i - 1];
        char opWhole = bam_cigar_op(valWhole);
        uint32_t lenWhole = bam_cigar_oplen(valWhole);

        if(opPart != opWhole) 
            return false;

        if(lenPart > lenWhole)
            return false;

        if(lenPart < lenWhole) {
            // length mismatch is only allowed in the last bases
            if(i != cigarNumPart-1) {
                // we only allow one CLIP in the end
                if(i != cigarNumPart-2)
                    return false;
                // check for next, is it CLIP ?
                int next = i+1;
                uint32_t valPartNext = cigarPart[next];
                // if it's right aligned, backward
                if(!isLeft)
                    valPartNext = cigarPart[cigarNumPart - next - 1];
                char opPartNext = bam_cigar_op(valPartNext);
                uint32_t lenPartNext = bam_cigar_oplen(valPartNext);
                if(opPartNext != BAM_CHARD_CLIP)
                    return false;
            }
        }
    }

    return true;
}

void BamUtil::dumpHeader(bam_hdr_t* hdr) {
    cerr << hdr->n_targets << " contigs in the bam file:" << endl;
    int dumped = 0;
    while(dumped < hdr->n_targets) {
        char *targetName = hdr->target_name[dumped];
        int targetLen = hdr->target_len[dumped];
        string name(targetName);
        cerr << targetName << ": " << targetLen << " bp" << endl;
        dumped++;
    }
    cerr << endl;
}

/* bam_cigar_type returns a bit flag with:
 *   bit 1 set if the cigar operation consumes the query
 *   bit 2 set if the cigar operation consumes the reference
 *
 * For reference, the unobfuscated truth table for this function is:
 * BAM_CIGAR_TYPE  QUERY  REFERENCE
 * --------------------------------
 * BAM_CMATCH      1      1
 * BAM_CINS        1      0
 * BAM_CDEL        0      1
 * BAM_CREF_SKIP   0      1
 * BAM_CSOFT_CLIP  1      0
 * BAM_CHARD_CLIP  0      0
 * BAM_CPAD        0      0
 * BAM_CEQUAL      1      1
 * BAM_CDIFF       1      1
 * BAM_CBACK       0      0
 * --------------------------------
 */

const int     QUERY_CONSUM[16] = {1, 1, 0, 0, 1, 0, 0, 1, 1, 0};
const int REFERENCE_CONSUM[16] = {1, 0, 1, 1, 0, 0, 0, 1, 1, 0};

int BamUtil::getRefOffset(bam1_t *b, int bampos) {
    uint32_t *data = (uint32_t *)bam_get_cigar(b);
    int cigarNum = b->core.n_cigar;
    int ref = 0;
    int query = 0;
    for(int i=0; i<cigarNum; i++) {
        uint32_t val = data[i];
        char op = bam_cigar_op(val);
        uint32_t len = bam_cigar_oplen(val);
        query += len * QUERY_CONSUM[op];
        ref += len * REFERENCE_CONSUM[op];
        if(query > bampos) {
            if(op == BAM_CINS || op == BAM_CSOFT_CLIP)
                return -1;
            else
                return ref - REFERENCE_CONSUM[op] * (query - bampos);
        }
    }
    cerr << "wrong cigar: " << getCigar(b) << " tid: " << b->core.tid << " l_qseq: " << b->core.l_qseq << " pos: "  << b->core.l_qseq << " isize: " << b->core.isize << " offset: " << bampos << endl;
    // not found
    return -1;
}

void BamUtil::getMOffsetAndLen(bam1_t *b, int& MOffset, int& MLen) {
    uint32_t *data = (uint32_t *)bam_get_cigar(b);
    int cigarNum = b->core.n_cigar;
    int ref = 0;
    int query = 0;
    for(int i=0; i<cigarNum; i++) {
        uint32_t val = data[i];
        char op = bam_cigar_op(val);
        uint32_t len = bam_cigar_oplen(val);
        if(op == BAM_CMATCH) {
            MOffset = query;
            MLen = len;
            return ;
        }
        query += len * QUERY_CONSUM[op];
        ref += len * REFERENCE_CONSUM[op];
    }
    // not found
    MOffset = 0;
    MLen = 0;
}

void BamUtil::copyQName(bam1_t *from, bam1_t *to) {
    char* fromdata = bam_get_qname(from);
    char* todata = bam_get_qname(to);
    int fromlen = from->core.l_qname;
    int tolen = to->core.l_qname;
    if(tolen < fromlen) {
        cerr << "copyQName ERROR: desitination qname is shorter";
        exit(-1);
    }

    memcpy(todata, fromdata, fromlen);

    // pad with \0
    for(int i=0; i< tolen - fromlen; i++) {
        todata[i + fromlen] = '\0';
    }
}

bool BamUtil::isPrimary(bam1_t *b) {
    if(b->core.flag & BAM_FSECONDARY || b->core.flag & BAM_FSUPPLEMENTARY)
        return false;
    else
        return true;
}

bool BamUtil::isProperPair(bam1_t *b) {
    return b->core.flag & BAM_FPROPER_PAIR;
}

int BamUtil::getRightRefPos(bam1_t *b) {
    if(b->core.pos<0)
        return -1;
    return b->core.pos + bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
}

bool BamUtil::test() {
    vector<string> qnames;
    qnames.push_back("NB551106:8:H5Y57BGX2:1:13304:3538:1404");
    qnames.push_back("NB551106:8:H5Y57BGX2:1:13304:3538:1404:UMI_GAGCATAC");
    qnames.push_back("NB551106:8:H5Y57BGX2:1:13304:3538:1404:UMI_GAGC_ATAC");
    qnames.push_back("NB551106:8:H5Y57BGX2:1:13304:3538:1404:GAGC_ATAC");
    qnames.push_back("NB551106:8:H5Y57BGX2:1:13304:3538:1404:UMI_X");

    vector<string> prefixes;
    prefixes.push_back("");
    prefixes.push_back("UMI");
    prefixes.push_back("UMI");
    prefixes.push_back("");
    prefixes.push_back("UMI");

    vector<string> umis;
    umis.push_back("");
    umis.push_back("GAGCATAC");
    umis.push_back("GAGC_ATAC");
    umis.push_back("GAGC_ATAC");
    umis.push_back("");

    for(int i=0; i<qnames.size(); i++) {
        string umi = getUMI(qnames[i], prefixes[i]);
        if(umi != umis[i]) {
            cerr << "get UMI from " << qnames[i] << ", expect " << umis[i] << ", but got " << umi << endl;
            return false;
        }
    }

    return true;

}
