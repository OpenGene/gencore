#include "pair.h"
#include "bamutil.h"

Pair::Pair(Options* opt){
    mLeft = NULL;
    mRight = NULL;
    mMergeReads = 1;
    mMergeLeftDiff = 0;
    mMergeRightDiff = 0;
    mOptions = opt;
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
}

void Pair::setLeft(bam1_t *b) {
    mLeft = b;
    mUMI = BamUtil::getUMI(mLeft, mOptions->umiPrefix);
    mLeftCigar = BamUtil::getCigar(mLeft);
}

void Pair::setRight(bam1_t *b) {
    mRight = b;
    mUMI = BamUtil::getUMI(mRight, mOptions->umiPrefix);
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
    cerr << "merged by " << mMergeReads << " reads, diff (" << mMergeLeftDiff << ", " << mMergeRightDiff << ")" << endl;
    if(mLeft)
        BamUtil::dump(mLeft);
    if(mRight)
        BamUtil::dump(mRight);
}

