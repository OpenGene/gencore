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

    cmdline::parser cmd;
    // input/output
    cmd.add<string>("in", 'i', "input sorted bam/sam file. STDIN will be read from if it's not specified", false, "-");
    cmd.add<string>("out", 'o', "output bam/sam file. STDOUT will be written to if it's not specified", false, "-");
    cmd.add<string>("ref", 'r', "reference fasta file name (should be an uncompressed .fa/.fasta file)", true, "");
    cmd.add<string>("umi_prefix", 'u', "the prefix for UMI, if it has. None by default. Check the README for the defails of UMI formats.", false, "");
    cmd.add<int>("supporting_reads", 's', "only output consensus reads/pairs that merged by >= <supporting_reads> reads/pairs. The valud should be 1~10, and the default value is 2.", false, 2);
    cmd.add<double>("ratio_threshold", 'a', "if the ratio of the major base in a cluster is less than <ratio_threshold>, it will be further compared to the reference. The valud should be 0.5~1.0, and the default value is 0.8", false, 0.8);
    cmd.add<int>("score_threshold", 'c', "if the score of the major base in a cluster is less than <score_threshold>, it will be further compared to the reference. The valud should be 1~20, and the default value is 6", false, 6);
    cmd.add<int>("quit_after_contig", 0, "stop when <quit_after_contig> contigs are processed. Only used for fast debugging. Default 0 means no limitation.", false, 0);
    cmd.add("debug", 0, "output some debug information to STDERR.");

    cmd.parse_check(argc, argv);

    Options opt;
    opt.input = cmd.get<string>("in");
    opt.output = cmd.get<string>("out");
    opt.refFile = cmd.get<string>("ref");
    opt.umiPrefix = cmd.get<string>("umi_prefix");
    opt.clusterSizeReq = cmd.get<int>("supporting_reads");
    opt.baseScoreReq = cmd.get<int>("score_threshold");
    opt.scorePercentReq = cmd.get<double>("ratio_threshold");
    opt.maxContig = cmd.get<int>("quit_after_contig");
    opt.debug = cmd.exist("debug");

    opt.validate();
    
    time_t t1 = time(NULL);

    // loading reference
    Reference* reference = NULL;
    if(!opt.refFile.empty()) {
        cerr << "loading reference data:" << endl;
        reference = Reference::instance(&opt);
    }

    Gencore gencore(&opt);
    gencore.consensus();

    if(reference) {
        delete reference;
        reference=NULL;
    }

    stringstream ss;
    for(int i=0;i<argc;i++){
        ss << argv[i] << " ";
    }
    command = ss.str();

    time_t t2 = time(NULL);
    cerr << endl << command << endl;
    cerr << "gencore v" << VERSION_NUMBER << ", time used: " << (t2)-t1 << " seconds" << endl;

    return 0;
}