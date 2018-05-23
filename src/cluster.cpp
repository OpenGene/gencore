#include "cluster.h"
#include "bamutil.h"
#include "reference.h"
#include <memory.h>

Cluster::Cluster(Options* opt){
    mOptions = opt;
}

Cluster::~Cluster(){
    for(int i=0; i<mPairs.size(); i++) {
        delete mPairs[i];
        mPairs[i] = NULL;
    }
}

void Cluster::addPair(Pair* p){
    mPairs.push_back(p);
}

void Cluster::dump(){
    for(int i=0; i<mPairs.size(); i++) {
        mPairs[i]->dump();
    }
}

bool Cluster::matches(Pair* p){
    for(int i=0; i<mPairs.size(); i++) {
        if(mPairs[i]->isDupWith(p))
        	return true;
    }
    return false;
}
    
vector<Pair*> Cluster::clusterByUMI(int umiDiffThreshold) {
	vector<Cluster*> subClusters;
    map<string, int> umiCount;
    for(int i=0; i<mPairs.size(); i++) {
        string umi = mPairs[i]->getUMI();
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

        Cluster* c = new Cluster(mOptions);

		// create the group by the top UMI
        vector<Pair*>::iterator piter;
        for(piter = mPairs.begin(); piter!=mPairs.end();){
    		Pair* p = *piter;
            string umi = p->getUMI();
            if(hamming(umi, topUMI) <= umiDiffThreshold) {
                c->addPair(p);
                piter = mPairs.erase(piter);
                umiCount[umi] = 0;
            } else {
                piter++;
            }
        }
        //if(mPairs.size()>0 || subClusters.size()>0)
        //    cerr << "UMI " << topUMI<< " " << topCount << "/" << c->mPairs.size() << endl;
        subClusters.push_back(c);
        umiCount[topUMI] = 0;
	}

    //if(subClusters.size()>1)
    //    cerr << subClusters.size() << " clusters" << endl;

	vector<Pair*> consensusPairs;

	for(int i=0; i<subClusters.size(); i++) {
		Pair* p = subClusters[i]->consensusMerge();
		consensusPairs.push_back(p);
		delete subClusters[i];
		subClusters[i] = NULL;
	}

	return consensusPairs;
}

Pair* Cluster::consensusMerge() {
    if(mPairs.size() == 1) {
    	Pair* p = mPairs[mPairs.size()-1];
    	p->mMergeReads = mPairs.size();
    	mPairs.pop_back();
    	return p;
    }

    int leftDiff = 0;
    int rightDiff = 0;
    bam1_t* left = consensusMergeBam(true, leftDiff);
    bam1_t* right = consensusMergeBam(false, rightDiff);

    Pair *p = new Pair(mOptions);
    p->mMergeReads = mPairs.size();
    if(left && right) {
        p->setLeft(left);
        p->mMergeLeftDiff = leftDiff;
        p->setRight(right);
        p->mMergeRightDiff = rightDiff;
        // make qname consistent to keep pair relationship
        if(left->core.l_qname <= right->core.l_qname) {
            BamUtil::copyQName(left, right);
        } else {
            BamUtil::copyQName(right, left);
        }
    }
    return p;
}

bam1_t* Cluster::consensusMergeBam(bool isLeft, int& diff) {
    bool leftReadMode = isLeft;
    // if processing right reads, check if this group is aligned by left
    if(!isLeft) {
        bool leftAligned = true;
        int lastPos = -1;
        for(int i=0; i<mPairs.size(); i++) {
            if(mPairs[i]->mRight) {
                if(lastPos >= 0 && mPairs[i]->mRight->core.pos != lastPos) {
                    leftAligned = false;
                    break;
                }
                lastPos = mPairs[i]->mRight->core.pos;
            }
        }
        // if it's left aligned, then process them as left reads
        if(leftAligned)
            leftReadMode = true;
    }
    // first we get a read that is most contained by other reads
    vector<int> containedByList(mPairs.size(), 0);
    for(int i=0; i<mPairs.size(); i++) {
        bam1_t* part = NULL;
        if(isLeft)
            part = mPairs[i]->mLeft;
        else
            part = mPairs[i]->mRight;
        if(part == NULL)
            continue;

        int containedBy = 0;

        for(int j=0; j<mPairs.size(); j++) {
            if(i == j)
                continue;
            bam1_t* whole = NULL;
            if(isLeft)
                whole = mPairs[j]->mLeft;
            else
                whole = mPairs[j]->mRight;
            if(whole == NULL)
                continue;

            // if processing right reads, we should align by the right ref pos
            if(!isLeft) {
                if(BamUtil::getRightRefPos(part) != BamUtil::getRightRefPos(whole)) {
                    continue;
                }
            }

            if( BamUtil::isPartOf(part, whole, leftReadMode))
                containedBy++;
        }

        containedByList[i] = containedBy;
    }

    int mostContainedById = -1;
    int mostContainedByNum = -1;
    for(int i=0; i<containedByList.size(); i++) {
        if(containedByList[i] > mostContainedByNum) {
            mostContainedByNum = containedByList[i];
            mostContainedById = i;
        } else if (containedByList[i] == mostContainedByNum && mostContainedById>=0) {
            // if they are identical, we use the shorter one
            int thisLen = 0;
            int curLen = 0;
            if(isLeft) {
                if(mPairs[i]->mLeft)
                    thisLen = mPairs[i]->mLeft->core.l_qseq;
                if(mPairs[mostContainedById]->mLeft)
                    curLen = mPairs[mostContainedById]->mLeft->core.l_qseq;
            } else {
                if(mPairs[i]->mRight)
                    thisLen = mPairs[i]->mRight->core.l_qseq;
                if(mPairs[mostContainedById]->mRight)
                    curLen = mPairs[mostContainedById]->mRight->core.l_qseq;
            }
            if(thisLen < curLen) {
                mostContainedByNum = containedByList[i];
                mostContainedById = i;
            }
        }
    }

    // no marjority
    if(mostContainedByNum < containedByList.size()*0.4)
        return NULL;

    bam1_t* out = NULL;
    if(isLeft) {
        out = mPairs[mostContainedById]->mLeft;
        // make it null so that it will not be deleted
        mPairs[mostContainedById]->mLeft = NULL;
    }
    else {
        out = mPairs[mostContainedById]->mRight;
        // make it null so that it will not be deleted
        mPairs[mostContainedById]->mRight = NULL;
    }

    vector<bam1_t *> reads;

    reads.push_back(out);

    for(int j=0; j<mPairs.size(); j++) {
        if(mostContainedById == j)
            continue;
        bam1_t* read = NULL;
        if(isLeft)
            read = mPairs[j]->mLeft;
        else
            read = mPairs[j]->mRight;
        if(read == NULL)
            continue;

        if( BamUtil::isPartOf(out, read, leftReadMode))
            reads.push_back(read);
    }

    if(reads.size() < mOptions->clusterSizeReq) {
        bam_destroy1(out);
        out = NULL;
        return NULL;
    }

    // if the sequences are right ones of pairs, we check whether it a completely a chaos
    /*if(!isLeft) {
        int bothEndNotAligned = 0;
        for(int r=0; r<reads.size(); r++) {
            // left aligned
            if(reads[r]->core.pos == out->core.pos && BamUtil::isPartOf(out, reads[r], true))
                continue;
            // right aligned
            if(reads[r]->core.pos + bam_cigar2rlen(reads[r]->core.n_cigar, (uint32_t *)bam_get_cigar(reads[r])) == out->core.pos + bam_cigar2rlen(out->core.n_cigar, (uint32_t *)bam_get_cigar(out)))
                continue;

            // both not aligned
            bothEndNotAligned++;
        }

        if(bothEndNotAligned*2 >= reads.size()) {
            cerr << "Chaos of " << reads.size() << " reads: " << BamUtil::getQName(out) << endl;
            for(int r=0; r<reads.size(); r++) {
                cerr << reads[r]->core.pos << "," << reads[r]->core.mpos << "," << reads[r]->core.isize << "," << reads[r]->core.l_qseq << "," << BamUtil::getCigar(reads[r]) << endl;
            }
            bam_destroy1(out);
            out = NULL;
            return NULL;
        }
    }*/

    diff = makeConsensus(reads, out, leftReadMode);

    return out;
}

int Cluster::makeConsensus(vector<bam1_t* >& reads, bam1_t* out, bool isLeft) {
    if(reads.size() <= 1 || out == NULL)
        return 0;

    int diff = 0;
    int mismatchInc = 0;

    // to restore the data if it's needed
    int seqbytes = (out->core.l_qseq+1)>>1;
    int qualbytes = out->core.l_qseq;
    char* seqBak = new char[seqbytes];
    char* qualBak = new char[qualbytes];
    memcpy(seqBak, bam_get_seq(out), seqbytes);
    memcpy(qualBak, bam_get_qual(out), qualbytes);

    vector<uint8_t *> alldata;
    vector<uint8_t *> allqual;
    vector<int> lenDiff;
    // if the sequences are right ones of pairs, we supposed they are aligned on the right (end)
    for(int r=0; r<reads.size(); r++) {
        alldata.push_back(bam_get_seq(reads[r]));
        allqual.push_back(bam_get_qual(reads[r]));
        int diff = reads[r]->core.l_qseq - out->core.l_qseq;
        if(diff != 0) {
            // A WAR for aligner induced unalignment on the right
            if(reads[r]->core.pos == out->core.pos && BamUtil::isPartOf(out, reads[r], true))
                diff = 0;
        }
        lenDiff.push_back(diff);
    }

    uint8_t * outdata = bam_get_seq(out);
    uint8_t * outqual = bam_get_qual(out);

    int len = out->core.l_qseq;
    if(out->core.n_cigar == 0) {
        for(int r=0; r<reads.size(); r++) {
            if(reads[r]->core.l_qseq < len)
                len = reads[r]->core.l_qseq;
        }
    }

    const char* refdata = NULL;
    if(out->core.isize != 0) {
        refdata = Reference::instance(mOptions)->getData(out->core.tid, out->core.pos, BamUtil::getRefOffset(out, len-1) + 1);
        if(refdata == NULL && mOptions->debug)
            cerr << "ref data is NULL for " << out->core.tid << ":" << out->core.pos << endl;
    }
    // loop all the position of out
    for(int i=0; i<len; i++) {
        int counts[16]={0};
        int quals[16]={0};
        int totalqual = 0;
        uint8_t topQual = 0;
        for(int r=0; r<reads.size(); r++) {
            int readpos = i;
            if(!isLeft)
                readpos = i + lenDiff[r];
            uint8_t base = 0;
            uint8_t qual = allqual[r][readpos];
            if(readpos%2 == 1)
                base = alldata[r][readpos/2] & 0xF;
            else
                base = (alldata[r][readpos/2]>>4) & 0xF;
            counts[base]++;
            quals[base] += qual;
            totalqual += qual;
            if(qual > topQual)
                topQual = qual;
        }
        // get the best representive base at this position
        uint8_t topBase=0;
        int topNum = 0;
        for(uint8_t b=0; b<16; b++) {
            if(counts[b] > topNum || (counts[b] == topNum && quals[b] > quals[topBase])) {
                topNum = counts[b];
                topBase = b;
            }
        }

        // get the secondary representive base at this position
        uint8_t secBase=0;
        int secNum = 0;
        for(uint8_t b=0; b<16; b++) {
            if(b == topBase)
                continue;
            if(counts[b] > secNum || (counts[b] == secNum && quals[b] > quals[secBase])) {
                secNum = counts[b];
                secBase = b;
            }
        }

        bool needToCheckRef = false;

        // completely consistent with good quality
        if(secNum == 0 && topQual >= 20) {
            outqual[i] = topQual;
            continue;
        }

        char refbase = 0;
        if(refdata) {
            int refpos = BamUtil::getRefOffset(out, i);
            if(refpos >= 0)
                refbase = refdata[refpos];
        }

        // the secondary base is a single base
        if(secNum ==1){
            // low quality secondary
            if(quals[secBase] < 15) {
                // the candidate is less than 3 consistent bases and no high quality bases
                if(topNum < 2 && topQual < 30) {
                    needToCheckRef = true;
                }
            }
            // high quality secondary
            else {
                // the candidate is less than 3 consistent bases or no high quality bases
                if(topNum < 3 || topQual < 30) {
                    needToCheckRef = true;
                }
            }
        }

        // more than one secondary bases
        if(secNum >1) {
            // the candidate is less than 80% consistent bases or no high quality
            if ((double)topNum < 0.8 * reads.size() || topQual < 30)
                needToCheckRef = true;
        }

        // integrate reference if it's possible
        if(needToCheckRef && refbase!=0) {
            uint8_t refbase4bit = BamUtil::base2fourbits(refbase);
            // check if there is one high quality base consistent to ref
            char refBaseQual = 0;
            for(int r=0; r<reads.size(); r++) {
                int readpos = i;
                if(!isLeft)
                    readpos = i + lenDiff[r];
                uint8_t base = 0;
                uint8_t qual = allqual[r][readpos];
                if(readpos%2 == 1)
                    base = alldata[r][readpos/2] & 0xF;
                else
                    base = (alldata[r][readpos/2]>>4) & 0xF;
                // found a ref-consistent base
                if(base == refbase4bit) {
                    // record the highest quality score of ref-consistent base
                    if(qual > refBaseQual)
                        refBaseQual = qual;
                    // if quality is high, use it for output
                    if(qual >= 30)
                        topBase = refbase4bit;
                }
            }
            // if there is no Q20 quality alternative, just use the reference base to reduce noise
            if(topQual < 20)
                topBase = refbase4bit;
            // if there is actually no ref base, the refBaseQual will have 0 quality
            if(topBase == refbase4bit)
                topQual = refBaseQual;
        }

        uint8_t outBase = 0;
        if(i%2 == 1)
            outBase = outdata[i/2] & 0xF;
        else
            outBase = (outdata[i/2]>>4) & 0xF;

        if(outBase != topBase) {
            // write this base to out
            if(i%2 == 1)
                outdata[i/2] = (outdata[i/2] & 0xF0) | topBase;
            else
                outdata[i/2] = (outdata[i/2] & 0x0F) | (topBase << 4);
            diff++;

            if(refbase!=0) {
                uint8_t refbase4bit = BamUtil::base2fourbits(refbase);
                if(outBase == refbase4bit) 
                    mismatchInc++;
                else if(topBase == refbase4bit) 
                    mismatchInc--;
            }
        }
        outqual[i] = topQual;
    }

    if(mismatchInc != 0) {

        // get edit distance info
        const char tagNM[2] ={'N', 'M'};
        uint8_t* dataNM = (uint8_t*)bam_aux_get(out,tagNM);
        int valNM = bam_aux2i(dataNM);
        int newValNM = valNM + mismatchInc;
        char typeNM = *dataNM;

        // this is abnormal, restore it, and output it for debugging
        if(mismatchInc > 5) {
            if(mOptions->debug) {
                cerr << endl;
                cerr << "NOTICE: mismatch increased with " << mismatchInc << endl;
                if(isLeft)
                    cerr << "Consensus by left" << endl;
                else
                    cerr << "Consensus by right" << endl;
                cerr << "Edit distance (NM) changed from " << valNM << " to " << newValNM << endl;
                cerr << "Read name: " << BamUtil::getQName(out) << endl;
                cerr << "tid: " << out->core.tid << ", pos: " << out->core.pos << endl;
                if(refdata)
                    cerr << "ref:" << string(refdata, len) << endl;
                cerr << "css:" << BamUtil::getSeq(out) << endl;
            }

            // restore the data
            memcpy(bam_get_seq(out),  seqBak,  seqbytes);
            memcpy(bam_get_qual(out), qualBak, qualbytes);

            if(mOptions->debug) {
                for(int r=0; r<reads.size(); r++) {
                    cerr << reads[r]->core.tid << ":" << reads[r]->core.pos << ", " << reads[r]->core.mpos << ", " << reads[r]->core.isize  << " " << BamUtil::getCigar(reads[r])  << " " << BamUtil::getSeq(reads[r]) << endl;
                }
                cerr << endl;
            }
        } else {
            // update edit distance
            if(typeNM == 'C' && newValNM>=0 && newValNM<=255) {
                dataNM[1] = newValNM;
            }
        }
    }

    delete seqBak;
    delete qualBak;

    return diff;
}

void Cluster::addRead(bam1_t* b) {
    // left
    if(b->core.isize > 0) {
        Pair* p = new Pair(mOptions);
        p->setLeft(b);
        mPairs.push_back(p);
        return;
    }
    // right or unproper paired
    bool found = false;
    string qname = BamUtil::getQName(b);
    for(int i=0; i<mPairs.size(); i++) {
        Pair* p = mPairs[i];
        if(p->mLeft && !p->mRight) {
            if(p->getQName() == qname) {
                p->setRight(b);
                found = true;
                break;
            }
        } else if(p->mRight && !p->mLeft) {
            if(p->getQName() == qname) {
                p->setLeft(b);
                found = true;
                break;
            }
        }
    }

    if(!found) {
        Pair* p = new Pair(mOptions);
        if(b->core.isize < 0)
            p->setRight(b);
        else
            p->setLeft(b);
        mPairs.push_back(p);
    }
}