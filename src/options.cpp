#include "options.h"
#include "util.h"

Options::Options(){
    input = "";
    output = "";
    refFile = "";
    umiPrefix = "";
    clusterSizeReq = 2;
    maxContig = 0;
    bamHeader = NULL;
    properReadsUmiDiffThreshold = 2;
    unproperReadsUmiDiffThreshold = 0;
    debug = false;

    scoreOfHighQualityMatch = 3;
    scoreOfLowQualityMatch = 2;
    scoreOfBothHighQualityMismatch = 0;
    scoreOfBothModerateQualityMismatch = -1;
    scoreOfBothLowQualityMismatch = -2;
    scoreOfUnbalancedMismatchHighQuality = 1;
    scoreOfUnbalancedMismatchLowQuality = -3;
}

bool Options::validate() {
    if(input.empty()) {
        error_exit("input should be specified by --in1");
    } else {
        check_file_valid(input);
    }

    if(ends_with(refFile, ".gz") || ends_with(refFile, ".gz")) {
        cerr << "reference fasta file should not be compressed.\nplease unzip "<<refFile<<" and try again."<<endl;
        exit(-1);
    }

    return true;
}