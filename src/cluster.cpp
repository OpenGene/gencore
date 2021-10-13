#include "cluster.h"
#include "bamutil.h"
#include "reference.h"
#include "group.h"
#include <memory.h>

Cluster::Cluster(Options* opt){
    mOptions = opt;
}

Cluster::~Cluster(){
    map<string, Pair*>::iterator iter;
    for(iter = mPairs.begin(); iter!=mPairs.end(); iter++) {
        delete iter->second;
    }
}

void Cluster::addPair(Pair* p){
    string qname = p->getQName();
    if(mPairs.count(qname)>0)
        delete mPairs[qname];
    mPairs[qname] = p;
}

void Cluster::dump(){
    map<string, Pair*>::iterator iter;
    for(iter = mPairs.begin(); iter!=mPairs.end(); iter++) {
        iter->second->dump();
    }
}

bool Cluster::matches(Pair* p){
    map<string, Pair*>::iterator iter;
    for(iter = mPairs.begin(); iter!=mPairs.end(); iter++) {
        if(iter->second->isDupWith(p))
            return true;
    }
    return false;
}

int Cluster::umiDiff(const string& umi1, const string& umi2) {

    int len1 = umi1.length();
    int len2 = umi2.length();

    int  diff = abs(len1 - len2);
    for(int i=0; i<min(len1, len2); i++) {
        if(umi1[i] != umi2[i])
             diff++;
    }

    return diff;
}
    
vector<Pair*> Cluster::clusterByUMI(int umiDiffThreshold, Stats* preStats, Stats* postStats, bool crossContig) {
	vector<Group*> groups;
    map<string, int> umiCount;
    bool hasUMI = false;
    map<string, Pair*>::iterator iterOfPairs;
    for(iterOfPairs = mPairs.begin(); iterOfPairs!=mPairs.end(); iterOfPairs++) {
        string umi = iterOfPairs->second->getUMI();
        if(!umi.empty())
            hasUMI = true;
        umiCount[umi]++;
    }
	while(mPairs.size()>0) {
        // get top UMI
        string topUMI;
        int topCount = 0;
        map<string, int>::iterator iter;
        for(iter = umiCount.begin(); iter!=umiCount.end(); iter++) {
            if(iter->second > topCount) {
                topCount = iter->second;
                topUMI = iter->first;
            }
        }

        Group* c = new Group(mOptions);
        bool isPE = false;

		// create the group by the top UMI
        map<string, Pair*>::iterator piter;
        for(piter = mPairs.begin(); piter!=mPairs.end();){
    		Pair* p = piter->second;
            string umi = p->getUMI();
            if(umiDiff(umi, topUMI) <= umiDiffThreshold) {
                c->addPair(p);
                if(p->mLeft && p->mRight)
                    isPE = true;
                piter = mPairs.erase(piter);
                umiCount[umi] = 0;
            } else {
                piter++;
            }
        }
        //if(mPairs.size()>0 || groups.size()>0)
        //    cerr << "UMI " << topUMI<< " " << topCount << "/" << c->mPairs.size() << endl;
        groups.push_back(c);
        umiCount[topUMI] = 0;
	}

    preStats->addCluster(groups.size()>1);

    //if(groups.size()>1)
    //    cerr << groups.size() << " clusters" << endl;

	vector<Pair*> singleConsensusPairs;

	for(int i=0; i<groups.size(); i++) {
		Pair* p = groups[i]->consensusMerge(crossContig);
		singleConsensusPairs.push_back(p);
		delete groups[i];
		groups[i] = NULL;
	}

    vector<Pair*> resultConsensusPairs;
    int singleConsesusCount = 0;
    int duplexConsensusCount = 0;
    if(hasUMI && !mOptions->disableDuplex) {
        // make duplex consensus read pairs
        while(singleConsensusPairs.size() > 0) {
            Pair* p1 = singleConsensusPairs.back();
            singleConsensusPairs.pop_back();
            string umi1  = p1->getUMI();
            bool foundDuplex = false;
            for(int i=0;i<singleConsensusPairs.size(); i++) {
                string umi2 = singleConsensusPairs[i]->getUMI();
                // a duplex
                if( isDuplex(umi1, umi2)) {
                    //cerr << "duplex:" << umi1 << ", " << umi2;
                    foundDuplex = true;
                    Pair* p2 = singleConsensusPairs[i];
                    // merge p2 to p1
                    int diff =  duplexMerge(p1, p2);
                    //cerr << " diff " << diff << endl;
                    preStats->addMolecule(p1->mMergeReads + p2->mMergeReads, p1->mLeft && p1->mRight);
                    if(diff <= mOptions->duplexMismatchThreshold) {
                        if(p1->mMergeReads + p2->mMergeReads >= mOptions->clusterSizeReq) {
                            duplexConsensusCount++;
                            p1->setDuplex(p2->mMergeReads);
                            p1->writeSscsDcsTag();
                            postStats->addDCS();
                            resultConsensusPairs.push_back(p1);
                        } else {
                            delete p1;
                        }
                    } else {
                        // too much diff, drop this duplex
                        delete p1;
                    }
                    singleConsensusPairs.erase(singleConsensusPairs.begin()+i);
                    delete p2;
                    break;
                }
            }
            // no duplex found, treat it as sscs
            if(!foundDuplex) {
                preStats->addMolecule(p1->mMergeReads, p1->mLeft && p1->mRight);
                if(!mOptions->duplexOnly && p1->mMergeReads >= mOptions->clusterSizeReq) {
                    singleConsesusCount++;
                    p1->writeSscsDcsTag();
                    postStats->addSSCS();
                    resultConsensusPairs.push_back(p1);
                } else {
                    delete p1;
                }
            }
        }
    } else {
        // no umi, no duplex
        for(int i=0;i<singleConsensusPairs.size(); i++) {
            Pair* p = singleConsensusPairs[i];
            preStats->addMolecule(p->mMergeReads, p->mLeft && p->mRight);
            if(!mOptions->duplexOnly && p->mMergeReads >= mOptions->clusterSizeReq) {
                singleConsesusCount++;
                p->writeSscsDcsTag();
                postStats->addSSCS();
                resultConsensusPairs.push_back(p);
            } else {
                delete p;
            }
        }
    }
    if(resultConsensusPairs.size()>0) {
        postStats->addCluster(resultConsensusPairs.size()>1);
    }
    return resultConsensusPairs;
}

int Cluster::duplexMerge(Pair* p1, Pair* p2) {
    int diff = 0;
    if(p1->mLeft && p2->mLeft)
        diff += duplexMergeBam(p1->mLeft, p2->mLeft);
    if(p1->mRight && p2->mRight)
        diff += duplexMergeBam(p1->mRight, p2->mRight);
    return diff;
}

int Cluster::duplexMergeBam(bam1_t* b1, bam1_t* b2) {
    int len1 = b1->core.l_qseq;
    int len2 = b2->core.l_qseq;
    int diff = abs(len1 - len2);
    int len = min(len1, len2);
    uint8_t * seq1 = bam_get_seq(b1);
    uint8_t * seq2 = bam_get_seq(b2);
    uint8_t * qual1 = bam_get_qual(b1);
    uint8_t * qual2 = bam_get_qual(b2);

    uint8_t N4bits = BamUtil::base2fourbits('N');
    for(int i=0; i<len; i++) {
        // two bases encoded in one byte: identical
        if(seq1[i/2] == seq2[i/2]) {
            i++;
            continue;
        }
        char base1, base2;
        if(i%2 == 1) {
            base1 = BamUtil::fourbits2base(seq1[i/2] & 0xF);
            base2 = BamUtil::fourbits2base(seq2[i/2] & 0xF);
        } else {
            base1 = BamUtil::fourbits2base((seq1[i/2]>>4) & 0xF);
            base2 = BamUtil::fourbits2base((seq2[i/2]>>4) & 0xF);
        }
        if(base1 != base2) {
            diff++;
            //set qual of the two bases to 0
            qual1[i]=0;
            qual2[i]=0;
            // set bases to N
            if(i%2 == 1) {
                seq1[i/2] &= 0xF0;
                seq1[i/2] |= N4bits;
                seq2[i/2] &= 0xF0;
                seq2[i/2] |= N4bits;
            } else {
                seq1[i/2] &= 0x0F;
                seq1[i/2] |= (N4bits << 4);
                seq2[i/2] &= 0x0F;
                seq2[i/2] |= (N4bits << 4);
            }
        }
    }
    return diff;
}

bool Cluster::isDuplex(const string& umi1, const string& umi2) {
    vector<string> umiPairs1;
    vector<string> umiPairs2;
    split(umi1, umiPairs1, "_");
    split(umi2, umiPairs2, "_");
    if(umiPairs1.size()!= 2 || umiPairs2.size() != 2)
        return false;

    if(umiPairs1[0] == umiPairs2[1] && umiPairs1[1] == umiPairs2[0])
        return true;
    else
        return false;
}

void Cluster::addRead(bam1_t* b) {
    // left
    string qname = BamUtil::getQName(b);
    map<string, Pair*>::iterator iter = mPairs.find(qname);

    if(iter!=mPairs.end()) {
        iter->second->setRight(b);
    }
    else {
        Pair* p = new Pair(mOptions);
        p->setLeft(b);
        mPairs[qname]=p;
    }
}

bool Cluster::test(){
    bool passed = true;
    passed &= umiDiff("ATCGATCG", "ATCGATCG") == 0;
    passed &= umiDiff("ATCGATCG", "ATCGTTC") == 2;
    passed &= umiDiff("ATCGATCG", "ATCGTTCG") == 1;
    passed &= umiDiff("AAAA_ATCG", "AAAA_ATCG") == 0;
    passed &= isDuplex("ATCG_CTAG", "CTAG_ATCG") == true;
    passed &= isDuplex("AGC_TGA", "TGA_AGC") == true;
    passed &= isDuplex("AAAA_AAAA", "AAAA_AAAA") == true;
    passed &= isDuplex("CTAG", "CTAG_ATCG") == false;
    passed &= isDuplex("CTAG", "CCCAGG") == false;
    passed &= isDuplex("", "") == false;
    return passed;
}