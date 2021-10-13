#include "pair.h"
#include "bamutil.h"
#include <memory.h>

Pair::Pair(Options* opt){
    mLeft = NULL;
    mRight = NULL;
    mLeftScore = NULL;
    mRightScore = NULL;
    mMergeReads = 1;
    mReverseMergeReads = 0;
    mMergeLeftDiff = 0;
    mMergeRightDiff = 0;
    mOptions = opt;
    mIsDuplex = false;
    mCssDcsTagWritten = false;
}

Pair::~Pair(){
    if(mLeft) {
        bam_destroy1(mLeft);
        mLeft = NULL;
    }
    if(mRight) {
        bam_destroy1(mRight);
        mRight = NULL;
    }
    if(mLeftScore) {
        delete[] mLeftScore;
        mLeftScore = NULL;
    }
    if(mRightScore) {
        delete[] mRightScore;
        mRightScore = NULL;
    }
}

void Pair::setDuplex(int mergeReadsOfReverseStrand) {
    mIsDuplex = true;
    mReverseMergeReads = mergeReadsOfReverseStrand;
}

void Pair::writeSscsDcsTag() {
    if(mCssDcsTagWritten) {
        error_exit("The SSCS/DCS tag has already been written!");
    }
    if(mLeft)
        writeSscsDcsTagBam(mLeft);
    if(mRight)
        writeSscsDcsTagBam(mRight);
    mCssDcsTagWritten = true;
}

void Pair::writeSscsDcsTagBam(bam1_t* b) {
    const char cssTag[2] = {'F','R'}; // forward strand read count
    const char dcsTag[2] = {'R','R'}; // reverse strand read number
    char type = 'C';
    unsigned short val = min(mMergeReads, 65535);
    int ret = bam_aux_append(b, cssTag, type, 1, (uint8_t*)&val);
    if(ret != 0)
        error_exit("Failed to write the consensus reads tag (CR) to BAM");
    if(mIsDuplex) {
        unsigned short valReverse = min(mReverseMergeReads, 65535);
        ret = bam_aux_append(b, dcsTag, type, 1, (uint8_t*)&valReverse);
        if(ret != 0)
            error_exit("Failed to write the duplex consensus reads tag (DR) to BAM");
    }
}

void Pair::assignNonOverlappedScores(uint8_t* qual, int start, int end, char* scores) {
    for(int i=start;i<end;i++) {
        uint8_t q = qual[i];
        scores[i] = qual2score(q);
    }
}

char Pair::qual2score(uint8_t q) {
    if(mOptions->highQuality <= q)
        return mOptions->scoreOfNotOverlappedHighQual;
    else if(mOptions->moderateQuality <= q)
        return mOptions->scoreOfNotOverlappedModerateQual;
    else if(mOptions->lowQuality <= q)
        return mOptions->scoreOfNotOverlappedLowQual;
    else
        return mOptions->scoreOfNotOverlappedBadQual;
}

void Pair::computeScore() {
    if(mLeft) {
        if(mLeftScore == NULL) {
            mLeftScore = new char[mLeft->core.l_qseq];
            memset(mLeftScore, mOptions->scoreOfNotOverlappedModerateQual, mLeft->core.l_qseq);
        }
    }

    if(mRight) {
        if(mRightScore == NULL) {
            mRightScore = new char[mRight->core.l_qseq];
            memset(mRightScore, mOptions->scoreOfNotOverlappedModerateQual, mRight->core.l_qseq);
        }
    }

    if(mLeftScore && mRightScore) {
        int leftMOffset, leftMLen, rightMOffset, rightMLen;
        BamUtil::getMOffsetAndLen(mLeft, leftMOffset, leftMLen);
        BamUtil::getMOffsetAndLen(mRight, rightMOffset, rightMLen);
        if(leftMLen>0 && rightMLen>0) {
            int posDis = mRight->core.pos - mLeft->core.pos;

            int leftStart, rightStart, cmpLen;
            if(posDis >=0 ) {
                leftStart = leftMOffset + posDis;
                rightStart = rightMOffset;
                cmpLen = min(leftMLen - posDis, rightMLen);
            } else {
                leftStart = leftMOffset;
                rightStart = rightMOffset - posDis;
                cmpLen = min(leftMLen, rightMLen + posDis);
            }
            uint8_t* lseq = bam_get_seq(mLeft);
            uint8_t* rseq = bam_get_seq(mRight);
            uint8_t* lqual = bam_get_qual(mLeft);
            uint8_t* rqual = bam_get_qual(mRight);
            if(mLeft) {
                assignNonOverlappedScores(lqual, 0, min(mLeft->core.l_qseq, leftStart), mLeftScore);
                assignNonOverlappedScores(lqual, max(0, leftStart+cmpLen), mLeft->core.l_qseq, mLeftScore);
            }
            if(mRight) {
                assignNonOverlappedScores(rqual, 0, min(mRight->core.l_qseq, rightStart), mRightScore);
                assignNonOverlappedScores(rqual, max(0, rightStart+cmpLen), mRight->core.l_qseq, mRightScore);
            }
            for(int i=0; i<cmpLen; i++) {
                int l = leftStart + i;
                int r = rightStart + i;
                uint8_t lbase, rbase;
                uint8_t lq = lqual[l];
                uint8_t rq = rqual[r];

                if(l%2 == 1)
                    lbase = lseq[l/2] & 0xF;
                else
                    lbase = (lseq[l/2]>>4) & 0xF;
                if(r%2 == 1)
                    rbase = rseq[r/2] & 0xF;
                else
                    rbase = (rseq[r/2]>>4) & 0xF;

                // matched pair, score += 4
                if(lbase == rbase) {
                    uint8_t q = (lq + rq) / 2;
                    char score = qual2score(q) + 4;
                    mLeftScore[l] = score;
                    mRightScore[r] = score;
                    continue;
                } else {
                    // modify the Q Scores since the bases are mismatched
                    // In the overlapped region, if a base and its pair are mismatched, its quality score will be adjusted to: max(0, this_qual - pair_qual)
                    lqual[l] = max(0, (int)lq - (int)rq);
                    rqual[r] = max(0, (int)rq - (int)lq);
                    uint8_t q = 0;
                    if(lq >= rq ) {
                        mLeftScore[l] = qual2score(lq - rq) - 3;
                        mRightScore[r] = 0;
                    } else {
                        mLeftScore[l] = 0;
                        mRightScore[r] = qual2score(rq - lq) -3;
                    }
                }
            }
        }
    }
}

char* Pair::getLeftScore() {
    if(!mLeftScore)
        computeScore();
    
    return mLeftScore;
}

char* Pair::getRightScore() {
    if(!mRightScore)
        computeScore();

    return mRightScore;
}

void Pair::setLeft(bam1_t *b) {
    if(mLeft)
        bam_destroy1(mLeft);
    mLeft = b;
    mUMI = BamUtil::getUMI(mLeft, mOptions->umiPrefix);
    mLeftCigar = BamUtil::getCigar(mLeft);
}

void Pair::setRight(bam1_t *b) {
    if(mRight)
        bam_destroy1(mRight);
    mRight = b;
    string umi = BamUtil::getUMI(mRight, mOptions->umiPrefix);
    if(!mUMI.empty() && umi!=mUMI) {
        cerr << "Mismatched UMI of a pair of reads" << endl;
        if(mLeft) {
            cerr << "Left:" << endl;
            BamUtil::dump(mLeft);
        }
        if(mRight) {
            cerr << "Right:" << endl;
            BamUtil::dump(mRight);
        }
        error_exit("The UMI of a read pair should be identical, but we got " + mUMI + " and " + umi );
    }
    else
        mUMI = umi;
    mRightCigar = BamUtil::getCigar(mRight);
}

bool Pair::pairFound() {
    return mLeft != NULL && mRight != NULL;
}

int Pair::getLeftRef() {
    if(mLeft == NULL)
        return -1;

    return mLeft->core.tid;
}

int Pair::getLeftPos() {
    if(mLeft == NULL)
        return -1;

    return mLeft->core.pos;
}

int Pair::getRightRef() {
    if(mRight == NULL)
        return -1;

    return mRight->core.tid;
}

int Pair::getRightPos() {
    if(mRight == NULL)
        return -1;

    return mRight->core.pos;
}

int Pair::getTLEN() {
    if(mLeft != NULL)
        return abs(mLeft->core.isize);
    else if(mRight != NULL)
        return abs(mRight->core.isize);
    else
        return -1;
}

string Pair::getQName() {
    if(mLeft != NULL)
        return BamUtil::getQName(mLeft);
    else if(mRight != NULL)
        return BamUtil::getQName(mRight);
    else
        return "";
}

Pair::MapType Pair::getMapType() {
    if(mLeft == NULL || mRight == NULL)
        return Unknown;

    int lref = getLeftRef();
    int rref = getRightRef();

    if(lref == rref) {
        if(lref >= 0)
            return ProperlyMapped;
        else
            return NoneMapped;
    }

    if(lref<0 && rref>=0)
        return OnlyRightMapped;

    if(lref>=0 && rref<0)
        return OnlyLeftMapped;

    if(lref != rref)
        return CrossRefMapped;

    return Unknown;
}

string Pair::getLeftCigar() {
    return mLeftCigar;
}

string Pair::getRightCigar() {
    return mRightCigar;
}

string Pair::getUMI() {
    return mUMI;
}


bool Pair::isDupWith(Pair* other) {
    if(!pairFound() || !other->pairFound())
        return false;

    if(getLeftRef() != other->getLeftRef() || getRightRef() != other->getRightRef())
        return false;

    if(getLeftPos() != other->getLeftPos() || getRightPos() != other->getRightPos())
        return false;

    int umiDiff = 0;
    for(int i=0; i<mUMI.length() && i<other->mUMI.length(); i++) {
        if(mUMI[i] != other->mUMI[i])
            umiDiff++;
    }
    if(umiDiff>1)
        return false;

    return true;
}

void Pair::dump() {
    cerr << "merged by " << mMergeReads << " forward reads, diff (" << mMergeLeftDiff << ", " << mMergeRightDiff << ")" << endl;
    if(mLeft){
        cerr << "left:" << endl;
        BamUtil::dump(mLeft);
    }
    if(mRight) {
        cerr << "right:" << endl;
        BamUtil::dump(mRight);
    }
}

