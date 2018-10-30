#ifndef JSON_REPORTER_H
#define JSON_REPORTER_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "options.h"
#include "stats.h"
#include <fstream>

using namespace std;

class JsonReporter{
public:
    JsonReporter(Options* opt);
    ~JsonReporter();

    void report(Stats* preStats, Stats* postStats);

private:
    Options* mOptions;
};


#endif