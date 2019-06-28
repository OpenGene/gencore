#include "bed.h"
#include "util.h"
#include <sstream>
#include <string.h>

Bed::Bed(Options* opt) {
	mOptions = opt;

	// initialize it
	// should be sure that bamHeader is already initialized
	for(int t=0; t < mOptions->bamHeader->n_targets;t++) {
        mContigRegions.push_back(vector<BedRegion>());
    }
}

void Bed::dump() {
	for(int c=0; c<mContigRegions.size(); c++) {
		for(int p=0; p<mContigRegions[c].size(); p++) {
			mContigRegions[c][p].dump();
		}
	}
}

vector<vector<long>> Bed::getDepthList() {
	vector<vector<long>> list;
	for(int c=0; c<mContigRegions.size(); c++) {
		list.push_back(vector<long>());
		for(int p=0; p<mContigRegions[c].size(); p++) {
			list[c].push_back(mContigRegions[c][p].getAvgDepth());
		}
	}
	return list;
}

string Bed::getPlotX(int c) {
	if(c >= mContigRegions.size())
		return "";

	stringstream ss;
	for(int p=0; p<mContigRegions[c].size(); p++) {
		ss << "\"" << mContigRegions[c][p].mName << " " << mContigRegions[c][p].mStart << "-" << mContigRegions[c][p].mEnd << "\"";
		if(p!= mContigRegions[c].size()-1)
			ss << ",";
	}
	return ss.str();
}

string Bed::getPlotY(int c, bool negative) {
	if(c >= mContigRegions.size())
		return "";

	stringstream ss;
	for(int p=0; p<mContigRegions[c].size(); p++) {
		if(negative)
			ss << "\"" << -mContigRegions[c][p].getAvgDepth() << "\"";
		else
			ss << "\"" << mContigRegions[c][p].getAvgDepth() << "\"";
		if(p!= mContigRegions[c].size()-1)
			ss << ",";
	}
	return ss.str();
}

void Bed::statDepth(int tid, int start, int len) {
	if(tid >= mContigRegions.size() || tid<0)
		return;

	int end = start + len;

	for(int p=0; p<mContigRegions[tid].size(); p++) {
		if(mContigRegions[tid][p].mEnd < start)
			continue;
		if(mContigRegions[tid][p].mStart > end)
			break;

		int len = min(mContigRegions[tid][p].mEnd, end) - max(mContigRegions[tid][p].mStart, start);
		mContigRegions[tid][p].mCount += len;
	}
}

void Bed::reportJSON(ofstream& ofs) {
	ofs << "," << endl;
	ofs << "\t\t\"coverage_bed\":{" << endl;
	for(int c=0; c<mContigRegions.size();c++) {
		string contig(mOptions->bamHeader->target_name[c]);
		ofs << "\t\t\t\"" << contig << "\":[" << endl;
		for(int p=0; p<mContigRegions[c].size(); p++) {
			ofs << "\t\t\t\t[\"" << mContigRegions[c][p].mName << "\"," << mContigRegions[c][p].mStart << "," << mContigRegions[c][p].mEnd << "," << mContigRegions[c][p].getAvgDepth() << "]";
			if(p != mContigRegions[c].size()-1)
				ofs << ",";
			ofs << endl;
		}
	    ofs << "\t\t\t]";
        if(c!=mContigRegions.size() - 1)
			ofs << ",";
        ofs << endl;
	}
	ofs << "\t\t}" << endl;
}

void Bed::copyFrom(Bed* other) {
	mContigRegions.clear();
	for(int c=0; c<other->mContigRegions.size(); c++) {
		mContigRegions.push_back(vector<BedRegion>());
		for(int p=0; p<other->mContigRegions[c].size(); p++) {
			mContigRegions[c].push_back(other->mContigRegions[c][p]);
		}
	}
}

void Bed::loadFromFile() {
	if(mOptions->bedFile.empty())
		return;
	else
		check_file_valid(mOptions->bedFile);

    ifstream file;
    file.open(mOptions->bedFile.c_str(), ifstream::in);
    const int maxLine = 4096;
    char line[maxLine];
    string lastChr;
    int lastTid;
    while(file.getline(line, maxLine)){
        // trim \n, \r or \r\n in the tail
        int readed = strlen(line);
        if(readed >=2 ){
            if(line[readed-1] == '\n' || line[readed-1] == '\r'){
                line[readed-1] = '\0';
                if(line[readed-2] == '\r')
                    line[readed-2] = '\0';
            }
        }
        string linestr(line);
        linestr = trim(linestr);
        vector<string> splitted;
        split(linestr, splitted, "\t");
        // comment line
        if(starts_with(splitted[0], "#"))
            continue;
        // require chr, start, end
        if(splitted.size()<3)
            continue;

        string chr = trim(splitted[0]);
        int start = atoi(trim(splitted[1]).c_str());
        int end = atoi(trim(splitted[2]).c_str());
        string name = "";
        if(splitted.size() > 3)
        	name = trim(splitted[3]);

        int tid = -1;
        if(chr == lastChr)
        	tid = lastTid;
        else {
        	for(int t=0; t < mOptions->bamHeader->n_targets;t++) {
		        char *targetName = mOptions->bamHeader->target_name[t];
		        string bamchr(targetName);
		        if(bamchr == chr) {
		        	tid = t;
		        	lastChr = chr;
		        	lastTid = tid;
		        }
		    }
        }

        if(tid>=0 && tid < mContigRegions.size())
        	mContigRegions[tid].push_back(BedRegion(chr, start, end, name));
    }
    mOptions->hasBedFile = true;
}