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

void Pair::computeScore() {
    if(mLeft) {
        if(mLeftScore == NULL) {
            mLeftScore = new char[mLeft->core.l_qseq];
            memset(mLeftScore, mOptions->scoreOfNotOverlapped, mLeft->core.l_qseq);
        }
    }

    if(mRight) {
        if(mRightScore == NULL) {
            mRightScore = new char[mRight->core.l_qseq];
            memset(mRightScore, mOptions->scoreOfNotOverlapped, mRight->core.l_qseq);
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

                if(lbase == rbase) {
                    if(lq + rq >= mOptions->moderateQuality * 2) {
                        mLeftScore[l] = mOptions->scoreOfHighQualityMatch;
                        mRightScore[r] = mOptions->scoreOfHighQualityMatch;
                    } else {
                        mLeftScore[l] = mOptions->scoreOfLowQualityMatch;
                        mRightScore[r] = mOptions->scoreOfLowQualityMatch;
                    }
                    continue;
                } else {
                    // modify the Q Scores since the bases are mismatched
                    lqual[l] = max(0, (int)lq - (int)rq);
                    rqual[r] = max(0, (int)rq - (int)lq);
                    if(lq >= mOptions->highQuality && rq >= mOptions->highQuality) {
                        mLeftScore[l] = mOptions->scoreOfBothHighQualityMismatch;
                        mRightScore[r] = mOptions->scoreOfBothHighQualityMismatch;
                        continue;
                    }

                    if(lq <= mOptions->lowQuality && rq <= mOptions->lowQuality) {
                        mLeftScore[l] = mOptions->scoreOfBothLowQualityMismatch;
                        mRightScore[r] = mOptions->scoreOfBothLowQualityMismatch;
                        continue;
                    }

                    if(lq > mOptions->lowQuality && lq < mOptions->highQuality && rq > mOptions->lowQuality && rq < mOptions->highQuality) {
                        mLeftScore[l] = mOptions->scoreOfBothModerateQualityMismatch;
                        mRightScore[r] = mOptions->scoreOfBothModerateQualityMismatch;
                        continue;
                    }

                    if(lq >= mOptions->highQuality && rq <= mOptions->lowQuality) {
                        mLeftScore[l] = mOptions->scoreOfUnbalancedMismatchHighQuality;
                        mRightScore[r] = mOptions->scoreOfUnbalancedMismatchLowQuality;
                        continue;
                    }

                    if(lq <= mOptions->lowQuality && rq >= mOptions->highQuality) {
                        mLeftScore[l] = mOptions->scoreOfUnbalancedMismatchLowQuality;
                        mRightScore[r] = mOptions->scoreOfUnbalancedMismatchHighQuality;
                        continue;
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

