#include "cluster.h"
#include "bamutil.h"
#include "reference.h"
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
	vector<Cluster*> subClusters;
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

        Cluster* c = new Cluster(mOptions);
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
        //if(mPairs.size()>0 || subClusters.size()>0)
        //    cerr << "UMI " << topUMI<< " " << topCount << "/" << c->mPairs.size() << endl;
        subClusters.push_back(c);
        preStats->addMolecule(c->mPairs.size(), isPE);
        umiCount[topUMI] = 0;
	}

    preStats->addCluster(subClusters.size()>1);

    //if(subClusters.size()>1)
    //    cerr << subClusters.size() << " clusters" << endl;

	vector<Pair*> singleConsensusPairs;

	for(int i=0; i<subClusters.size(); i++) {
		Pair* p = subClusters[i]->consensusMerge(crossContig);
		singleConsensusPairs.push_back(p);
		delete subClusters[i];
		subClusters[i] = NULL;
	}

    vector<Pair*> resultConsensusPairs;
    int singleConsesusCount = 0;
    int duplexConsensusCount = 0;
    if(hasUMI) {
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
                if(p1->mMergeReads >= mOptions->clusterSizeReq) {
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
            if(p->mMergeReads >= mOptions->clusterSizeReq) {
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
            base2 = BamUtil::fourbits2base((seq1[i/2]>>4) & 0xF);
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

Pair* Cluster::consensusMerge(bool crossContig) {
    int leftDiff = 0;
    int rightDiff = 0;

    // in this case, no need to make consensus
    if(mPairs.size()==1 && mPairs.begin()->second->mRight == NULL) {
        Pair* p = mPairs.begin()->second;
        mPairs.clear();
        return p;
    }

    bam1_t* nameToCopy = NULL;
    if(crossContig) {
        int curLen = 0;

        map<string, Pair*>::iterator iter;
        for(iter=mPairs.begin(); iter!=mPairs.end(); iter++) {
            if(iter->second->mLeft == NULL)
                continue;

            if(nameToCopy == NULL) {
                nameToCopy = iter->second->mLeft;
                curLen = iter->second->mLeft->core.l_qname;
                continue;
            }

            if(iter->second->mLeft->core.l_qname < curLen || (iter->second->mLeft->core.l_qname == curLen && strcmp(bam_get_qname(iter->second->mLeft), bam_get_qname(nameToCopy)) <0 )) {
                nameToCopy = iter->second->mLeft;
                curLen = iter->second->mLeft->core.l_qname;
            }
        }
    }

    bam1_t* left = consensusMergeBam(true, leftDiff);
    bam1_t* right = consensusMergeBam(false, rightDiff);

    Pair *p = new Pair(mOptions);
    p->mMergeReads = mPairs.size();

    // for cross-contig mapped reads, only left read is present
    // to keep the PE relationship, we use the smallest name in the shortest read names
    if(crossContig) {
        if(left && nameToCopy && nameToCopy != left) {
            BamUtil::copyQName(nameToCopy, left);
        }

    } else if(left && right) {
        string lname = BamUtil::getQName(left);
        string rname = BamUtil::getQName(right);
        // make qname consistent to keep pair relationship
        if( lname.length() <=  rname.length()) {
            BamUtil::copyQName(left, right);
        } else {
            BamUtil::copyQName(right, left);
        }
    }
    if(left) {
        p->setLeft(left);
        p->mMergeLeftDiff = leftDiff;
    }

    if(right) {
        p->setRight(right);
        p->mMergeRightDiff = rightDiff;
    }
    return p;
}

bam1_t* Cluster::consensusMergeBam(bool isLeft, int& diff) {
    vector<Pair*> allPairs;
    map<string, Pair*>::iterator iterOfPairs;
    for(iterOfPairs = mPairs.begin(); iterOfPairs!=mPairs.end(); iterOfPairs++) {
        allPairs.push_back(iterOfPairs->second);
    }
    if(mPairs.size() > mOptions->skipLowComplexityClusterThreshold) {
        map<string, int> cigars;
        bam1_t* firstRead = NULL;
        for(iterOfPairs = mPairs.begin(); iterOfPairs!=mPairs.end(); iterOfPairs++) {
            Pair* p = iterOfPairs->second;
            bam1_t* b = p->mLeft;
            if(!isLeft)
                b = p->mRight;
            if(b) {
                string cigar = BamUtil::getCigar(b);
                if(cigars.count(cigar) == 0)
                    cigars[cigar] = 1;
                else
                    cigars[cigar]++;
                if(!firstRead)
                    firstRead = b;
            }
        }
        // this is abnormal, usually due to mapping result of low complexity reads
        if(cigars.size() > mPairs.size() * 0.1 && firstRead) {
            string seq = BamUtil::getSeq(firstRead);
            int diffNeighbor = 0;
            for(int i=0;i<seq.length()-1;i++) {
                if(seq[i] != seq[i+1])
                    diffNeighbor++;
            }
            if(diffNeighbor < seq.length()*0.5) {
                if(mOptions->debug) {
                    cerr << "Skipping " << mPairs.size() << " low complexity reads like: " << seq << endl;
                }
                return NULL;
            }
        }
    }

    bool leftReadMode = isLeft;
    // if processing right reads, check if this group is aligned by left
    if(!isLeft) {
        bool leftAligned = true;
        int lastPos = -1;
        for(int i=0; i<allPairs.size(); i++) {
            if(allPairs[i]->mRight) {
                if(lastPos >= 0 && allPairs[i]->mRight->core.pos != lastPos) {
                    leftAligned = false;
                    break;
                }
                lastPos = allPairs[i]->mRight->core.pos;
            }
        }
        // if it's left aligned, then process them as left reads
        if(leftAligned)
            leftReadMode = true;
    }
    // first we get a read that is most contained by other reads
    vector<int> containedByList(allPairs.size(), 0);
    for(int i=0; i<allPairs.size(); i++) {
        bam1_t* part = NULL;
        if(isLeft)
            part = allPairs[i]->mLeft;
        else
            part = allPairs[i]->mRight;
        if(part == NULL)
            continue;

        int containedBy = 1;

        for(int j=0; j<allPairs.size(); j++) {
            if(i == j)
                continue;
            bam1_t* whole = NULL;
            if(isLeft)
                whole = allPairs[j]->mLeft;
            else
                whole = allPairs[j]->mRight;
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
        if(mPairs.size() > mOptions->skipLowComplexityClusterThreshold && containedBy>=mPairs.size()/2) 
            break;
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
                if(allPairs[i]->mLeft)
                    thisLen = allPairs[i]->mLeft->core.l_qseq;
                if(allPairs[mostContainedById]->mLeft)
                    curLen = allPairs[mostContainedById]->mLeft->core.l_qseq;
            } else {
                if(allPairs[i]->mRight)
                    thisLen = allPairs[i]->mRight->core.l_qseq;
                if(allPairs[mostContainedById]->mRight)
                    curLen = allPairs[mostContainedById]->mRight->core.l_qseq;
            }
            if(thisLen < curLen) {
                mostContainedByNum = containedByList[i];
                mostContainedById = i;
            }
        }
    }

    // no marjority
    if(mostContainedByNum < containedByList.size()*0.4 && containedByList.size() != 1) {
        return NULL;
    }

    bam1_t* out = NULL;
    char* outScore = NULL;
    if(isLeft) {
        out = allPairs[mostContainedById]->mLeft;
        outScore = allPairs[mostContainedById]->getLeftScore();
        // make it null so that it will not be deleted
        allPairs[mostContainedById]->mLeft = NULL;
    }
    else {
        out = allPairs[mostContainedById]->mRight;
        outScore = allPairs[mostContainedById]->getRightScore();
        // make it null so that it will not be deleted
        allPairs[mostContainedById]->mRight = NULL;
    }

    if(out == NULL) {
        return NULL;
    }

    vector<bam1_t *> reads;
    vector<char *> scores;

    reads.push_back(out);
    scores.push_back(outScore);

    for(int j=0; j<allPairs.size(); j++) {
        if(mostContainedById == j)
            continue;
        bam1_t* read = NULL;
        char* score = NULL;
        if(isLeft) {
            read = allPairs[j]->mLeft;
            score = allPairs[j]->getLeftScore();
        }
        else {
            read = allPairs[j]->mRight;
            score = allPairs[j]->getRightScore();
        }
        if(read == NULL || score == NULL)
            continue;

        if( BamUtil::isPartOf(out, read, leftReadMode)) {
            reads.push_back(read);
            scores.push_back(score);
        }
    }

    diff = makeConsensus(reads, out, scores, leftReadMode);

    return out;
}

int Cluster::makeConsensus(vector<bam1_t* >& reads, bam1_t* out, vector<char*>& scores, bool isLeft) {
    if(out == NULL)
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

    const unsigned char* refdata = NULL;
    if(out->core.isize != 0) {
        refdata = Reference::instance(mOptions)->getData(out->core.tid, out->core.pos, BamUtil::getRefOffset(out, len-1) + 1);
        if(refdata == NULL && mOptions->debug)
            cerr << "ref data is NULL for " << out->core.tid << ":" << out->core.pos << endl;
    }
    // loop all the position of out
    for(int i=0; i<len; i++) {
        int counts[16]={0};
        int baseScores[16]={0};
        int quals[16]={0};
        uint8_t topQuals[16] = {0};
        int totalqual = 0;
        int totalScore = 0;
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
            baseScores[base] += scores[r][readpos];
            totalScore += scores[r][readpos];
            quals[base] += qual;
            totalqual += qual;
            if(qual > topQuals[base])
                topQuals[base] = qual;
        }
        // get the best representive base at this position
        uint8_t topBase=0;
        int topScore = -0x7FFFFFFF;
        for(uint8_t b=0; b<16; b++) {
            if(baseScores[b] > topScore || (baseScores[b] == topScore && quals[b] > quals[topBase])) {
                topScore = baseScores[b];
                topBase = b;
            }
        }
        int topNum = counts[topBase];
        uint8_t topQual = topQuals[topBase];

        // get the secondary representive base at this position
        uint8_t secBase=0;
        int secScore = -0x7FFFFFFF;
        for(uint8_t b=0; b<16; b++) {
            if(b == topBase)
                continue;
            if(baseScores[b] > secScore || (baseScores[b] == secScore && quals[b] > quals[secBase])) {
                secScore = baseScores[b];
                secBase = b;
            }
        }
        int secNum = counts[secBase];

        bool needToCheckRef = false;

        // if the secondary base is not a valid candidate and the top base is valid
        if(secNum == 0) {
            if(topScore >= mOptions->baseScoreReq && topQual >= mOptions->moderateQuality) {
                outqual[i] = topQual;
                continue;
            } else
                needToCheckRef = true;
        }

        char refbase = 0;
        if(refdata) {
            int refpos = BamUtil::getRefOffset(out, i);
            if(refpos >= 0) {
                refbase = FastaReader::getBase(refdata, out->core.pos + refpos);
            }
        }

        if(refbase!='A' && refbase!='T' && refbase!='C' && refbase!='G')
            refbase = 0;

        // the secondary base is a single base
        if(secNum ==1){
            // low quality secondary
            if(quals[secBase] <= mOptions->lowQuality) {
                // the candidate is less than 2 consistent bases and no high quality bases
                if(topNum < 2 && topQual < mOptions->highQuality) {
                    needToCheckRef = true;
                }
            }
            // high quality secondary
            else {
                // the candidate is less than 3 consistent bases or no high quality bases
                if(topNum < 3 || topQual < mOptions->highQuality) {
                    needToCheckRef = true;
                }
            }
        }

        // more than one secondary bases
        if(secNum >1) {
            // the candidate is less than 80% consistent bases or they are all bad quality
            if ((double)topScore < mOptions->scorePercentReq * totalScore || topQual < mOptions->moderateQuality)
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
                    if(qual >= mOptions->highQuality)
                        topBase = refbase4bit;
                }
            }
            // if there is no alternative with moderate quality, just use the reference base to reduce noise
            if(topQual < mOptions->moderateQuality)
                topBase = refbase4bit;
            // if there is actually no ref base, the refBaseQual will have 0 quality
            // so that it will be masked for downstream processing
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
                    cerr << "ref:" << endl << FastaReader::toString(refdata, out->core.pos, len) << endl;
                cerr << "css:" << endl << BamUtil::getSeq(out) << endl;
            }

            // restore the data
            memcpy(bam_get_seq(out),  seqBak,  seqbytes);
            memcpy(bam_get_qual(out), qualBak, qualbytes);

            if(mOptions->debug) {
                for(int r=0; r<reads.size(); r++) {
                    cerr << reads[r]->core.tid << ":" << reads[r]->core.pos << ", " << reads[r]->core.mpos << ", " << reads[r]->core.isize  << " " << BamUtil::getCigar(reads[r])  << endl << BamUtil::getSeq(reads[r]) << endl;
                    for(int p=0; p<reads[r]->core.l_qseq; p++)
                        cerr << (int)scores[r][p];
                    cerr<<endl;
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

    delete[] seqBak;
    delete[] qualBak;

    return diff;
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