#include "options.h"
#include "util.h"

Options::Options(){
    input = "";
    output = "";
    refFile = "";
    umiPrefix = "";
    maxContig = 0;
    bamHeader = NULL;
    properReadsUmiDiffThreshold = 2;
    unproperReadsUmiDiffThreshold = 0;
    debug = false;

    scorePercentReq = 0.8;
    clusterSizeReq = 2;
    baseScoreReq = 1;

    highQuality = 30;
    moderateQuality = 20;
    lowQuality = 15;

    scoreOfNotOverlapped = 1;
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

    if(scorePercentReq > 1.0) {
        error_exit("ratio_threshold cannot be greater than 1.0");
    } else if(scorePercentReq < 0.5) {
        error_exit("ratio_threshold cannot be less than 0.5");
    }

    if(clusterSizeReq > 10) {
        error_exit("supporting_reads cannot be greater than 10");
    } else if(clusterSizeReq < 1) {
        error_exit("supporting_reads cannot be less than 1");
    }

    if(baseScoreReq > 10) {
        error_exit("score_threshold cannot be greater than 10");
    } else if(baseScoreReq < 1) {
        error_exit("score_threshold cannot be less than 1");
    }

    return true;
}