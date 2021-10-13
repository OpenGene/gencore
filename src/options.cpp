#include "options.h"
#include "util.h"

Options::Options(){
    input = "";
    output = "";
    refFile = "";
    bedFile = "";
    umiPrefix = "";
    maxContig = 0;
    bamHeader = NULL;
    properReadsUmiDiffThreshold = 1;
    unproperReadsUmiDiffThreshold = 0;
    duplexMismatchThreshold = 2;
    debug = false;
    hasBedFile = false;

    scorePercentReq = 0.8;
    clusterSizeReq = 1;

    highQuality = 30;
    moderateQuality = 20;
    lowQuality = 15;

    scoreOfNotOverlappedHighQual = 8;
    scoreOfNotOverlappedModerateQual = 6;
    scoreOfNotOverlappedLowQual = 4;
    scoreOfNotOverlappedBadQual = 2;

    baseScoreReq = scoreOfNotOverlappedModerateQual;
    skipLowComplexityClusterThreshold = 1000;

    reportTitle = "gencore report";

    bedCoverageStep = 10;
    coverageStep = 10000;

    duplexOnly = false;
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

    if(highQuality > 40) {
        error_exit("high_qual cannot be greater than 40");
    } else if(highQuality < 20) {
        error_exit("high_qual cannot be less than 20");
    }

    if(moderateQuality > 35) {
        error_exit("moderate_qual cannot be greater than 35");
    } else if(moderateQuality < 15) {
        error_exit("moderate_qual cannot be less than 15");
    }

    if(lowQuality > 30) {
        error_exit("low_qual cannot be greater than 30");
    } else if(lowQuality < 8) {
        error_exit("low_qual cannot be less than 8");
    }

    if(properReadsUmiDiffThreshold > 10) {
        error_exit("umi_diff_threshold cannot be greater than 10");
    } else if(properReadsUmiDiffThreshold < 0) {
        error_exit("umi_diff_threshold cannot be negative");
    }

    if(lowQuality > moderateQuality) {
        error_exit("low_qual cannot be greater than moderate_qual");
    }

    if(moderateQuality > highQuality) {
        error_exit("moderate_qual cannot be greater than high_qual");
    }

    if(duplexMismatchThreshold > 10) {
        error_exit("duplex_diff_threshold cannot be greater than 10, suggest 2.");
    } else if(duplexMismatchThreshold < 0) {
        error_exit("duplex_diff_threshold cannot be less than 0, suggest 2.");
    }

    return true;
}