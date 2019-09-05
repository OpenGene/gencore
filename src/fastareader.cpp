
#include "fastareader.h"
#include "util.h"
#include <sstream>
#include <memory>
#include <string.h>

FastaReader::FastaReader(Options* opt, string faFile, bool forceUpperCase)
{
    // Set locale and disable stdio synchronization to improve iostream performance
    // http://www.drdobbs.com/the-standard-librarian-iostreams-and-std/184401305
    // http://stackoverflow.com/questions/5166263/how-to-get-iostream-to-perform-better
    mOptions = opt;
    setlocale(LC_ALL,"C");
    ios_base::sync_with_stdio(false);

    mFastaFile = faFile;
    mForceUpperCase = forceUpperCase;
    if (is_directory(mFastaFile)) {
        string error_msg = "There is a problem with the provided fasta file: \'";
        error_msg.append(mFastaFile);
        error_msg.append("\' is a directory NOT a file...\n");
        throw invalid_argument(error_msg);
    }
    mFastaFileStream.open( mFastaFile.c_str(),ios::in);
    // verify that the file can be read
    if (!mFastaFileStream.is_open()) {
        string msg = "There is a problem with the provided fasta file: could NOT read ";
        msg.append(mFastaFile.c_str());
        msg.append("...\n");
        throw invalid_argument(msg);
    }

    char c;
    // seek to first contig
    while (mFastaFileStream.get(c) && c != '>') {
        if (mFastaFileStream.eof()) {
            break;
        }
    }
}

FastaReader::~FastaReader()
{
    if (mFastaFileStream.is_open()) {
        mFastaFileStream.close();
    }

    map<string, unsigned char*>::iterator iter;
    for(iter = mAllContigs.begin(); iter != mAllContigs.end(); iter++) {
        if(iter->second) {
            delete[] iter->second;
            iter->second = NULL;
        }
    }
}

void FastaReader::readNext()
{
    mCurrentID = "";
    mCurrentDescription = "";
    mCurrentSequence = NULL;
    mCurrentSize = 0;
    bool foundHeader = false;
    
    char c;
    stringstream ssSeq;
    stringstream ssHeader;
    while(true){
        mFastaFileStream.get(c);
        if(c == '>' || mFastaFileStream.eof())
            break;
        else {
            if (foundHeader){
                if(mForceUpperCase && c>='a' && c<='z') {
                    c -= ('a' - 'A');
                }
                ssSeq << c;
            }
            else
                ssHeader << c;
        }

        string line = "";
        getline(mFastaFileStream,line,'\n');


        if(foundHeader == false) {
            ssHeader << line;
            foundHeader = true;
        }
        else {
            str_keep_valid_sequence(line, mForceUpperCase);
            ssSeq << line;
        }
    }
    string str = ssSeq.str();
    mCurrentSize = str.length();
    mCurrentSequence = to4bits(str);
    string header = ssHeader.str();

    int space = header.find(" ");
    mCurrentID = header.substr(0, space);
}

unsigned char FastaReader::base2bits(char base) {
    if(base =='A') return 1;
    else if(base =='T') return 2;
    else if(base =='C') return 3;
    else if(base =='G') return 4;
    else
        return 0; // N or others
}

char FastaReader::bits2base(unsigned char bits) {
    if(bits >= 5)
        return 'N';
    const char bases[5] ={'N', 'A', 'T', 'C', 'G'};
    return bases[bits];
}

unsigned char FastaReader::getBase(const unsigned char* refdata, int refpos) {
    char two4bits = refdata[refpos/2];
    if(refpos % 2 ==0)
        return bits2base(two4bits & 0x0F);
    else
        return bits2base( (two4bits & 0xF0) >> 4 );
}

string FastaReader::toString(const unsigned char* refdata, int pos, int len) {
    string str(len, 'N');
    for(int i=0; i<len; i++)
        str[i] = getBase(refdata, pos+i);

    return str;
}

// encode two four bits base in one byte
unsigned char* FastaReader::to4bits(const string & str) {
    size_t len = (str.length() + 1) / 2;
    unsigned char* data = new unsigned char[len];
    memset(data, 0, len);
    for(int i=0; i<str.length(); i++) {
        char val = 0;
        char bits = base2bits(str[i]);
        if(i%2 == 0)
            data[i/2] |= bits;
        else
            data[i/2] |= (bits << 4);
    }
    return data;
}

bool FastaReader::hasNext() {
    return !mFastaFileStream.eof();
}

void FastaReader::readAll() {
    while(!mFastaFileStream.eof()){
        readNext();
        cerr << mCurrentID << ": " << mCurrentSize << " bp" << endl;
        mAllContigs[mCurrentID] = mCurrentSequence;
        mAllContigSizes[mCurrentID] = mCurrentSize;
        if(mOptions->maxContig>0 && mAllContigs.size()>mOptions->maxContig){
            break;
        }
    }
    cerr << endl << "loaded " << mAllContigs.size() << " contigs" << endl<< endl;
}



