#include "stats.h"
#include <memory.h>
#include <sstream>
#include "util.h"


Stats::Stats(Options* opt) {
	mOptions = opt;
	memset(this, 0, sizeof(Stats));
	mSupportingHistgram = new long[MAX_SUPPORTING_READS];
	memset(mSupportingHistgram, 0, sizeof(long)*MAX_SUPPORTING_READS);
}

Stats::~Stats() {
	delete[] mSupportingHistgram;
}

long Stats::getMappedBases() {
	return mBase - mBaseUnmapped;
}

long Stats::getMappedReads() {
	return mRead - mReadUnmapped;
}

long Stats::getReads() {
	return mRead;
}

long Stats::getBases() {
	return mBase;
}

void Stats::addRead(int baseNum, int mismatch, bool mapped) {
	mBase += baseNum;
	mRead++;
	mBaseMismatches += mismatch;

	if(!mapped) {
		mBaseUnmapped += baseNum;
		mReadUnmapped++;
	}

	if(mismatch>0)
		mReadWithMismatches++;
}

void Stats::addMolecule(unsigned int supportingReads, bool PE) {
	mMolecule++;
	if(supportingReads < MAX_SUPPORTING_READS)
		mSupportingHistgram[supportingReads]++;
	if(PE)
		mMoleculePE++;
	else
		mMoleculeSE++;
}

void Stats::addCluster(bool hasMultiMolecule) {
	mCluster++;
	if(hasMultiMolecule)
		mMultiMoleculeCluster++;
}

double Stats::getMappingRate() {
	return getMappedReads() / (double)mRead;
}

double Stats::getDupRate() {
	return 1.0 - (mMoleculeSE + mMoleculePE * 2) / (double) getMappedReads();
}

double Stats::getMismatchRate() {
	return (double)mBaseMismatches/getMappedBases();
}

void Stats::reportJSON(ofstream& ofs) {
	ofs << "\t\t\"total_reads\": " << mRead << "," << endl;
	ofs << "\t\t\"total_bases\": " << mBase << "," << endl;
	ofs << "\t\t\"mapped_reads\": " << getMappedReads() << "," << endl;
	ofs << "\t\t\"mapped_bases\": " << getMappedBases() << "," << endl;
	ofs << "\t\t\"mismatched_bases\": " << mBaseMismatches  <<  "," << endl;
	ofs << "\t\t\"reads_with_mismatched_bases\": " << mReadWithMismatches <<  "," << endl;
	ofs << "\t\t\"mismatch_rate\": "  <<  (double)mBaseMismatches/getMappedBases() << "," << endl;
	ofs << "\t\t\"total_mapping_clusters\": " << mCluster << "," << endl;
	ofs << "\t\t\"multiple_fragments_clusters\": " << mMultiMoleculeCluster << "," << endl;
	ofs << "\t\t\"total_fragments\": " << mMolecule << "," << endl;
	ofs << "\t\t\"single_end_fragments\": " << mMoleculeSE << "," << endl;
	ofs << "\t\t\"paired_end_fragments\": " << mMoleculePE << "," << endl;
	ofs << "\t\t\"duplication_level_histogram\": [";
	for(int i=1; i<MAX_SUPPORTING_READS - 1; i++)
		ofs << mSupportingHistgram[i] << ",";
	ofs << mSupportingHistgram[MAX_SUPPORTING_READS-1];
	ofs << "]" << endl;
}

void Stats::print() {
	cerr << "Total reads: " << mRead << endl;
	cerr << "Total bases: " << mBase << endl;
	cerr << "Mapped reads: " << getMappedReads() << " (" << to_string( getMappedReads()*100.0/mRead ) << "%)" << endl;
	cerr << "Mapped bases: " << getMappedBases() << " (" << to_string( getMappedBases()*100.0/mBase ) << "%)" << endl;
	cerr << "Bases mismatched with reference: " << mBaseMismatches << " (" << to_string( mBaseMismatches*100.0/getMappedBases() ) << "%)" <<  endl;
	cerr << "Reads with mismatched bases: " << mReadWithMismatches << " (" << to_string( mReadWithMismatches*100.0/getMappedReads() ) << "%)" <<  endl;
	cerr << "Total mapping clusters: " << mCluster << endl;
	cerr << "Mapping clusters with multiple fragments: " << mMultiMoleculeCluster << endl;
	cerr << "Total fragments: " << mMolecule << endl;
	cerr << "Fragments with single-end reads: " << mMoleculeSE << endl;
	cerr << "Fragments with paired-end reads: " << mMoleculePE << endl;
	cerr << "Duplication level histogram: " << endl;
	for(int i=1; i<MAX_SUPPORTING_READS && i<=10; i++) {
		if(mSupportingHistgram[i] == 0)
			break;
		cerr << "    Fragments with " << i << " duplicates: " << mSupportingHistgram[i] << endl;
	}
}

string Stats::list2string(double* list, int size) {
    stringstream ss;
    for(int i=0; i<size; i++) {
        ss << list[i];
        if(i < size-1)
            ss << ",";
    }
    return ss.str();
}

string Stats::list2string(double* list, int size, long* coords) {
    stringstream ss;
    for(int i=0; i<size; i++) {
        // coords is 1,2,3,...
        long start = 0;
        if(i>0)
            start = coords[i-1];
        long end = coords[i];

        double total = 0.0;
        for(int k=start; k<end; k++)
            total += list[k];

        // get average
        if(end == start)
            ss << "0.0";
        else
            ss << total / (end - start);
        //ss << list[coords[i]-1];
        if(i < size-1)
            ss << ",";
    }
    return ss.str();
}

string Stats::list2string(long* list, int size) {
    stringstream ss;
    for(int i=0; i<size; i++) {
        ss << list[i];
        if(i < size-1)
            ss << ",";
    }
    return ss.str();
}