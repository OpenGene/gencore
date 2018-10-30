#include "jsonreporter.h"

JsonReporter::JsonReporter(Options* opt){
    mOptions = opt;
}

JsonReporter::~JsonReporter(){
}

extern string command;
void JsonReporter::report(Stats* preStats, Stats* postStats) {
	ofstream ofs;
    ofs.open(mOptions->jsonFile, ifstream::out);
    ofs << "{" << endl;

    // before
    ofs << "\t" << "\"summary\": {" << endl;
    ofs << "\t\t\"mapping_rate\":" << preStats->getMappingRate() << "," << endl;
    ofs << "\t\t\"duplication_rate\":" << preStats->getDupRate() << "," << endl;
    ofs << "\t\t\"fragment_passing_rate\":" << (double) postStats->getMolecules() / preStats->getMolecules() << "";
    ofs << endl;
    ofs << "\t" << "}," << endl;

    // before
    ofs << "\t" << "\"before_processing\": {" << endl;
    preStats->reportJSON(ofs);
    ofs << endl;
    ofs << "\t" << "}," << endl;

    // after
    ofs << "\t" << "\"after_processing\": {" << endl;
    postStats->reportJSON(ofs);
    ofs << endl;
    ofs << "\t" << "}," << endl;

    ofs << "\t\"command\": " << "\"" << command << "\"" << endl;

    ofs << "}";

    ofs.close();
}