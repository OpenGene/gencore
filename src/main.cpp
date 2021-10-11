#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "cmdline.h"
#include "common.h"
#include <sstream>
#include "util.h"
#include "gencore.h"
#include "options.h"
#include "reference.h"
#include "unittest.h"

using namespace std;

string command;

int main(int argc, char* argv[]){
    if (argc == 2 && strcmp(argv[1], "test")==0){
        UnitTest tester;
        tester.run();
        return 0;
    }

    if (argc == 2 && (strcmp(argv[1], "-v")==0 || strcmp(argv[1], "--version")==0)){
        cerr << "gencore " << VERSION_NUMBER << endl;
        return 0;
    }

    cmdline::parser cmd;
    // input/output
    cmd.add<string>("in", 'i', "input sorted bam/sam file. STDIN will be read from if it's not specified", false, "-");
    cmd.add<string>("out", 'o', "output bam/sam file. STDOUT will be written to if it's not specified", false, "-");
    cmd.add<string>("ref", 'r', "reference fasta file name (should be an uncompressed .fa/.fasta file)", true, "");
    cmd.add<string>("bed", 'b', "bed file to specify the capturing region, none by default", false, "");
    cmd.add("duplex_only", 'x', "only output duplex consensus sequences, which means single stranded consensus sequences will be discarded.");
    
    // UMI
    cmd.add<string>("umi_prefix", 'u', "the prefix for UMI, if it has. None by default. Check the README for the defails of UMI formats.", false, "auto");
    
    // thresholds
    cmd.add<int>("supporting_reads", 's', "only output consensus reads/pairs that merged by >= <supporting_reads> reads/pairs. The valud should be 1~10, and the default value is 1.", false, 1);
    cmd.add<double>("ratio_threshold", 'a', "if the ratio of the major base in a cluster is less than <ratio_threshold>, it will be further compared to the reference. The valud should be 0.5~1.0, and the default value is 0.8", false, 0.8);
    cmd.add<int>("score_threshold", 'c', "if the score of the major base in a cluster is less than <score_threshold>, it will be further compared to the reference. The valud should be 1~20, and the default value is 6", false, 6);
    cmd.add<int>("umi_diff_threshold", 'd', "if two reads with identical mapping position have UMI difference <= <umi_diff_threshold>, then they will be merged to generate a consensus read. Default value is 1.", false, 1);
    cmd.add<int>("duplex_diff_threshold", 'D', "if the forward consensus and reverse consensus sequences have <= <duplex_diff_threshold> mismatches, then they will be merged to generate a duplex consensus sequence, otherwise will be discarded. Default value is 2.", false, 2);
    cmd.add<int>("high_qual", 0, "the threshold for a quality score to be considered as high quality. Default 30 means Q30.", false, 30);
    cmd.add<int>("moderate_qual", 0, "the threshold for a quality score to be considered as moderate quality. Default 20 means Q20.", false, 20);
    cmd.add<int>("low_qual", 0, "the threshold for a quality score to be considered as low quality. Default 15 means Q15.", false, 15);
    cmd.add<int>("coverage_sampling", 0, "the sampling rate for genome scale coverage statistics. Default 10000 means 1/10000.", false, 10000);

    // reporting
    cmd.add<string>("json", 'j', "the json format report file name", false, "gencore.json");
    cmd.add<string>("html", 'h', "the html format report file name", false, "gencore.html");

    // debugging
    cmd.add("debug", 0, "output some debug information to STDERR.");
    cmd.add<int>("quit_after_contig", 0, "stop when <quit_after_contig> contigs are processed. Only used for fast debugging. Default 0 means no limitation.", false, 0);

    cmd.parse_check(argc, argv);

    Options opt;
    opt.input = cmd.get<string>("in");
    opt.output = cmd.get<string>("out");
    opt.refFile = cmd.get<string>("ref");
    opt.bedFile = cmd.get<string>("bed");
    opt.umiPrefix = cmd.get<string>("umi_prefix");
    opt.clusterSizeReq = cmd.get<int>("supporting_reads");
    opt.baseScoreReq = cmd.get<int>("score_threshold");
    opt.scorePercentReq = cmd.get<double>("ratio_threshold");
    opt.maxContig = cmd.get<int>("quit_after_contig");
    opt.highQuality = cmd.get<int>("high_qual");
    opt.moderateQuality = cmd.get<int>("moderate_qual");
    opt.lowQuality = cmd.get<int>("low_qual");
    opt.coverageStep = cmd.get<int>("coverage_sampling");
    opt.properReadsUmiDiffThreshold = cmd.get<int>("umi_diff_threshold");
    opt.duplexMismatchThreshold = cmd.get<int>("duplex_diff_threshold");
    opt.debug = cmd.exist("debug");
    opt.duplexOnly = cmd.exist("duplex_only");

    // reporting
    opt.jsonFile = cmd.get<string>("json");
    opt.htmlFile = cmd.get<string>("html");

    opt.validate();
    
    time_t t1 = time(NULL);

    // loading reference
    Reference* reference = NULL;
    if(!opt.refFile.empty()) {
        cerr << "loading reference data:" << endl;
        reference = Reference::instance(&opt);
    }

    stringstream ss;
    for(int i=0;i<argc;i++){
        ss << argv[i] << " ";
    }
    command = ss.str();

    Gencore gencore(&opt);
    gencore.consensus();

    if(reference) {
        delete reference;
        reference=NULL;
    }

    time_t t2 = time(NULL);
    cerr << endl << command << endl;
    cerr << "gencore v" << VERSION_NUMBER << ", time used: " << (t2)-t1 << " seconds" << endl;

    return 0;
}