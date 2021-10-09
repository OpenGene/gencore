#include "htmlreporter.h"
#include <chrono>
#include <memory.h>
#include <math.h>
#include "common.h"

extern string command;

HtmlReporter::HtmlReporter(Options* opt){
    mOptions = opt;
    mDupRate = 0.0;
}

HtmlReporter::~HtmlReporter(){
}

void HtmlReporter::setInsertHist(long* insertHist, int insertSizePeak) {
    mInsertHist = insertHist;
    mInsertSizePeak = insertSizePeak;
}

void HtmlReporter::outputRow(ofstream& ofs, string key, long v) {
    ofs << "<tr><td class='col1'>" + key + "</td><td class='col2'>" + to_string(v) + "</td></tr>\n";
}

void HtmlReporter::outputRow(ofstream& ofs, string key, string v) {
    ofs << "<tr><td class='col1'>" + key + "</td><td class='col2'>" + v + "</td></tr>\n";
}

void HtmlReporter::outputTripleRow(ofstream& ofs, string key, string v1, string v2) {
    ofs << "<tr><td class='col1'>" + key + "</td><td class='col2'>" + v1 + "</td><td class='col3'>" + v2 + "</td></tr>\n";
}

string HtmlReporter::formatNumber(long number) {
    double num = (double)number;
    string unit[6] = {"", "K", "M", "G", "T", "P"};
    int order = 0;
    while (num > 1000.0) {
        order += 1;
        num /= 1000.0;
    }

    if (order == 0)
        return to_string(number);
    else
        return to_string(num) + " " + unit[order];
}

string HtmlReporter::getPercents(long numerator, long denominator) {
    if(denominator == 0)
        return "0.0";
    else
        return to_string((double)numerator * 100.0 / (double)denominator);
}

void HtmlReporter::printSummary(ofstream& ofs,  Stats* preStats, Stats* postStats) {

    ofs << endl;
    ofs << "<h1 style='text-align:left;'><a href='https://github.com/OpenGene/gencore' target='_blank' style='color:#663355;text-decoration:none;'>" + mOptions->reportTitle + "</a>"<<endl;
    ofs << "<div class='section_div'>\n";
    ofs << "<div class='section_title' onclick=showOrHide('summary')><a name='summary'>Summary</a></div>\n";
    ofs << "<div id='summary'>\n";

    ofs << "<div class='subsection_title' onclick=showOrHide('general')>General</div>\n";
    ofs << "<div id='general'>\n";
    ofs << "<table class='summary_table'>\n";
    outputRow(ofs, "gencore version:", string(VERSION_NUMBER)+ " (<a href='https://github.com/OpenGene/gencore'>https://github.com/OpenGene/gencore</a>)");
    outputRow(ofs, "mapping rate:", to_string(preStats->getMappingRate()));
    outputRow(ofs, "duplication rate:", to_string(preStats->getDupRate()));
    outputRow(ofs, "Single Stranded Consensus Sequence:", to_string(postStats->mSSCSNum));
    outputRow(ofs, "Duplex Consensus Sequence:", to_string(postStats->mDCSNum));
    ofs << "</table>\n";
    ofs << "</div>\n";

    ofs << "<div class='subsection_title' onclick=showOrHide('filtering_metrics')>Details</div>\n";
    ofs << "<div id='filtering_metrics'>\n";
    ofs << "<table class='summary_table'>\n";
    outputTripleRow(ofs, "", "before processing", "after processing");
    outputTripleRow(ofs, "total bases:", formatNumber(preStats->getBases()),formatNumber(postStats->getBases()) );
    outputTripleRow(ofs, "mapped bases:", formatNumber(preStats->getMappedBases()),formatNumber(postStats->getMappedBases()) );
    outputTripleRow(ofs, "total reads:", formatNumber(preStats->getReads()),formatNumber(postStats->getReads()) );
    outputTripleRow(ofs, "mapped reads:", formatNumber(preStats->getMappedReads()),formatNumber(postStats->getMappedReads()) );
    outputTripleRow(ofs, "mismatched bases:", formatNumber(preStats->mBaseMismatches),formatNumber(postStats->mBaseMismatches) );
    outputTripleRow(ofs, "reads with mismatched bases:", formatNumber(preStats->mBaseMismatches),formatNumber(postStats->mBaseMismatches) );
    outputTripleRow(ofs, "mismatch rate:", to_string(preStats->getMismatchRate()),to_string(postStats->getMismatchRate()) );
    outputTripleRow(ofs, "total mapping clusters:", formatNumber(preStats->mCluster),formatNumber(postStats->mCluster) );
    outputTripleRow(ofs, "multiple fragments clusters:", formatNumber(preStats->mMultiMoleculeCluster),formatNumber(postStats->mMultiMoleculeCluster) );
    outputTripleRow(ofs, "total fragments:", formatNumber(preStats->mMolecule),formatNumber(postStats->mMolecule) );
    outputTripleRow(ofs, "single-end fragments:", formatNumber(preStats->mMoleculeSE),formatNumber(postStats->mMoleculeSE) );
    outputTripleRow(ofs, "paired-end fragments:", formatNumber(preStats->mMoleculePE),formatNumber(postStats->mMoleculePE) );
    ofs << "</table>\n";
    ofs << "</div>\n";

    ofs << "</div>\n";
    ofs << "</div>\n";

    if(true) {
        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('duplication')><a name='duplication'>Duplication histogram of mapped reads</a></div>\n";
        ofs << "<div id='duplication'>\n";

        reportDuplication(ofs, preStats->mSupportingHistgram, preStats);

        ofs << "</div>\n";
        ofs << "</div>\n";
    }

    if(true) {
        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('coverage')><a name='coverage'>Coverage statistics in genome scale</a></div>\n";
        ofs << "<div id='coverage'>\n";

        reportCoverage(ofs, preStats, postStats);

        ofs << "</div>\n";
        ofs << "</div>\n";
    }

    if(mOptions->hasBedFile) {
        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('coverage_bed')><a name='coverage_bed'>Coverage statistics in BED:<font size=-2>" << mOptions->bedFile << "</font> </a></div>\n";
        ofs << "<div id='coverage_bed'>\n";

        reportCoverageBed(ofs, preStats, postStats);

        ofs << "</div>\n";
        ofs << "</div>\n";
    }

    if(false) {
        ofs << "<div class='section_div'>\n";
        ofs << "<div class='section_title' onclick=showOrHide('insert_size')><a name='insert_size'>Insert size estimation</a></div>\n";
        ofs << "<div id='insert_size'>\n";

        //reportInsertSize(ofs, preStats1->getCycles() + preStats2->getCycles() - mOptions->overlapRequire);

        ofs << "</div>\n";
        ofs << "</div>\n";
    }
}

long HtmlReporter::getYCeiling(vector<vector<long>> list, int denominator) {
    int size = 0;
    for(int i=0; i<list.size(); i++) {
        if(mOptions->maxContig == 0 || i <= mOptions->maxContig)
            size += list[i].size();
    }
    size = 1 + size/denominator;
    long* topvalues = new long[size];
    memset(topvalues, 0, sizeof(long) * size);
    for(int x=0; x<list.size(); x++) {
        for(int y=0; y<list[x].size();y++) {
            long v = list[x][y];
            for(int j=size-1; j>=0; j--) {
                if(v > topvalues[j]) {
                    for(int p=0; p<j;p++) {
                        topvalues[p] = topvalues[p+1];
                    }
                    topvalues[j] = v;
                    break;
                }
            }
        }
    }
    return topvalues[0];
}

void HtmlReporter::reportCoverage(ofstream& ofs, Stats* preStats, Stats* postStats) {
    int maxpos = 0;
    for(int c=0; c<preStats->mGenomeDepth.size();c++) {
        if(preStats->mGenomeDepth[c].size() > maxpos) {
            maxpos = preStats->mGenomeDepth[c].size();
        }
    }
    // to avoid outliers in the figure
    double ceilingY = (double)getYCeiling(preStats->mGenomeDepth, 500)/mOptions->coverageStep;

    ofs << "<div style='padding:5px;'><center><table style='border:0px;'><tr><td style='width:20px;background:red'></td><td style='border:0px;'>Before processing</td><td style='width:20px;background:blue'></td><td style='border:0px;'>After processing</td></tr></table></center></div>"<<endl;

    for(int c=0; c<preStats->mGenomeDepth.size();c++) {

        if(preStats->mGenomeDepth[c].size() * 100 < maxpos)
            continue;

        double w = 5.0 + 95.0 * preStats->mGenomeDepth[c].size() / maxpos;
        string contig(mOptions->bamHeader->target_name[c]);

        double* ybefore = new double[preStats->mGenomeDepth[c].size()];
        double* yafter = new double[preStats->mGenomeDepth[c].size()];
        double* x = new double[preStats->mGenomeDepth[c].size()];

        int total = preStats->mGenomeDepth[c].size();
        for(int i=0; i< preStats->mGenomeDepth[c].size(); i++) {
            x[i] = i*mOptions->coverageStep;
            ybefore[i] = (double)preStats->mGenomeDepth[c][i] / mOptions->coverageStep;
            yafter[i] = -(double)postStats->mGenomeDepth[c][i] / mOptions->coverageStep;
        }


        ofs << "<div class='coverage_div' id='coverage_" + contig +"'>\n";
        ofs << "<div class='coverage_figure' id='plot_coverage_" + contig + "' style='width:" + to_string(w) + "%;height:80px;'></div>\n";
        ofs << "</div>\n";
        
        ofs << "\n<script type=\"text/javascript\">" << endl;
        string json_str = "var data=[";

        // original
        json_str += "{";
        json_str += "x:[" + Stats::list2string(x, total) + "],";
        json_str += "y:[" + Stats::list2string(ybefore, total) + "],";
        json_str += "name: 'before processing',";
        json_str += "fill: 'tozeroy',";
        json_str += "line:{color:'rgb(255,0, 0)', width:1}\n";
        json_str += "},";

        // original
        json_str += "{";
        json_str += "x:[" + Stats::list2string(x, total) + "],";
        json_str += "y:[" + Stats::list2string(yafter, total) + "],";
        json_str += "name: 'after processing',";
        json_str += "fill: 'tozeroy',";
        json_str += "line:{color:'rgb(0, 0, 255)', width:1}\n";
        json_str += "}";

        json_str += "];\n";

        json_str += "var layout={margin: {l: 50,r: 50,b:30,t: 5,pad: 2}, showlegend: false, yaxis:{title:'" + contig + "', range:[" + to_string(-ceilingY) + ", " + to_string(ceilingY)+ "]}};\n";
        json_str += "Plotly.newPlot('plot_coverage_" + contig + "', data, layout);\n";

        ofs << json_str;
        ofs << "</script>" << endl;

        delete[] ybefore;
        delete[] yafter;
        delete[] x;
    }
}

void HtmlReporter::reportCoverageBed(ofstream& ofs, Stats* preStats, Stats* postStats) {
    vector<vector<BedRegion>>& preBed = preStats->mBedStats->mContigRegions;
    vector<vector<BedRegion>>& postBed = postStats->mBedStats->mContigRegions;

    int maxpos = 0;
    for(int c=0; c<preBed.size();c++) {
        if(preBed[c].size() > maxpos) {
            maxpos = preBed[c].size();
        }
    }
    // to avoid outliers in the figure
    long ceilingY1 = getYCeiling(preStats->mBedStats->getDepthList(), 500);
    long ceilingY2 = getYCeiling(postStats->mBedStats->getDepthList(), 500);

    ofs << "<div style='padding:5px;'><center><table style='border:0px;'><tr><td style='width:20px;background:red'></td><td style='border:0px;'>Before processing</td><td style='width:20px;background:blue'></td><td style='border:0px;'>After processing</td></tr></table></center></div>"<<endl;

    for(int c=0; c<preBed.size();c++) {

        if(preBed[c].size() == 0)
            continue;

        double w = 5.0 + 95.0 * max(maxpos/100.0, (double)preBed[c].size()) / maxpos;
        string contig(mOptions->bamHeader->target_name[c]);

        double* ybefore = new double[preStats->mGenomeDepth[c].size()];
        double* yafter = new double[preStats->mGenomeDepth[c].size()];
        double* x = new double[preStats->mGenomeDepth[c].size()];

        int total = preBed[c].size();

        ofs << "<div class='bed_coverage_div' id='bed_coverage_" + contig +"'>\n";
        ofs << "<div class='coverage_figure' id='bed_plot_coverage_" + contig + "' style='width:" + to_string(w) + "%;height:250px;'></div>\n";
        ofs << "</div>\n";
        
        ofs << "\n<script type=\"text/javascript\">" << endl;
        string json_str = "var data=[";

        // original
        json_str += "{";
        json_str += "x:[" + preStats->mBedStats->getPlotX(c) + "],";
        json_str += "y:[" + preStats->mBedStats->getPlotY(c) + "],";
        json_str += "name: 'before processing',";
        json_str += "fill: 'tozeroy',";
        json_str += "line:{color:'rgb(255,0, 0)', width:1}\n";
        json_str += "},";

        // original
        json_str += "{";
        json_str += "x:[" + postStats->mBedStats->getPlotX(c) + "],";
        json_str += "y:[" + postStats->mBedStats->getPlotY(c, true) + "],";
        json_str += "name: 'after processing',";
        json_str += "fill: 'tozeroy',";
        json_str += "line:{color:'rgb(0, 0, 255)', width:1}\n";
        json_str += "}";

        json_str += "];\n";

        json_str += "var layout={margin: {l: 50,r: 50,b: 150,t:5,pad: 2}, xaxis:{tickangle:60, tickfont:{size: 8,color: '#bc6f98'}}, showlegend: false, yaxis:{title:'" + contig + "', range:[" + to_string(-ceilingY2) + ", " + to_string(ceilingY1)+ "]}};\n";
        json_str += "Plotly.newPlot('bed_plot_coverage_" + contig + "', data, layout);\n";

        ofs << json_str;
        ofs << "</script>" << endl;

    }
}

void HtmlReporter::reportInsertSize(ofstream& ofs, int isizeLimit) {
    if(isizeLimit<1)
        isizeLimit = 1;
    int insertSizeMax = 200;
    int total = min(insertSizeMax, isizeLimit);
    long *x = new long[total];
    double allCount = 0;
    for(int i=0; i<total; i++) {
        x[i] = i;
        allCount += mInsertHist[i];
    }
    allCount += mInsertHist[insertSizeMax];
    double* percents = new double[total];
    memset(percents, 0, sizeof(double)*total);
    if(allCount > 0) {
        for(int i=0; i<total; i++) {
            percents[i] = (double)mInsertHist[i] * 100.0 / (double)allCount;
        }
    }

    double unknownPercents = (double)mInsertHist[insertSizeMax] * 100.0 / (double)allCount;

    ofs << "<div id='insert_size_figure'>\n";
    ofs << "<div class='figure' id='plot_insert_size' style='height:400px;'></div>\n";
    ofs << "</div>\n";
    
    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";

    json_str += "{";
    json_str += "x:[" + Stats::list2string(x, total) + "],";
    json_str += "y:[" + Stats::list2string(percents, total) + "],";
    json_str += "name: 'Percent (%)  ',";
    json_str += "";
    json_str += "line:{color:'rgba(128,0,128,1.0)', width:1}\n";
    json_str += "}";

    json_str += "];\n";

    json_str += "var layout={title:'Insert size distribution (" + to_string(unknownPercents) + "% reads are with unknown length)', xaxis:{title:'Insert size'}, yaxis:{title:'Read percent (%)'}};\n";
    json_str += "Plotly.newPlot('plot_insert_size', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;

    delete[] x;
    delete[] percents;
}

void HtmlReporter::reportDuplication(ofstream& ofs, long* dupHist, Stats* preStats) {

    ofs << "<div id='duplication_figure'>\n";
    ofs << "<div class='figure' id='plot_duplication' style='height:400px;'></div>\n";
    ofs << "</div>\n";
    
    ofs << "\n<script type=\"text/javascript\">" << endl;
    string json_str = "var data=[";

    int total = MAX_SUPPORTING_READS - 1;
    while(total > 0 && dupHist[total] == 0) {
        total--;
    }
    if(total == 0)
        total = 1;
    long *x = new long[total];
    double allCount = preStats->uncountedSupportingReads;
    for(int i=0; i<total; i++) {
        x[i] = i+1;
        allCount += dupHist[i+1];
    }
    double* percents = new double[total];
    memset(percents, 0, sizeof(double)*total);
    double uncountedPercent = 0.0; 
    if(allCount > 0) {
        for(int i=0; i<total; i++) {
            percents[i] = (double)dupHist[i+1] * 100.0 / (double)allCount;
        }
        uncountedPercent = 100.0 * preStats->uncountedSupportingReads / allCount;
    }

    json_str += "{type:'bar',";
    json_str += "x:[" + Stats::list2string(x, total) + "],";
    json_str += "y:[" + Stats::list2string(percents, total) + "],";
    json_str += "name: 'Read percent (%)  ',";
    json_str += "";
    json_str += "line:{color:'rgba(128,0,128,1.0)'}\n";
    json_str += "},";

    json_str += "];\n";

    json_str += "var layout={title:'" + to_string(uncountedPercent) + " % fragments have " + to_string(MAX_SUPPORTING_READS) + "+ duplicated reads', xaxis:{title:'duplication level'}, yaxis:{title:'Fragment percent (%)'}};\n";
    json_str += "Plotly.newPlot('plot_duplication', data, layout);\n";

    ofs << json_str;
    ofs << "</script>" << endl;

    delete[] x;
    delete[] percents;
}

void HtmlReporter::report(Stats* preStats, Stats* postStats) {
    ofstream ofs;
    ofs.open(mOptions->htmlFile, ifstream::out);

    printHeader(ofs);

    printSummary(ofs,  preStats, postStats);

    printFooter(ofs);

}

void HtmlReporter::printHeader(ofstream& ofs){
    ofs << "<html><head><meta http-equiv=\"content-type\" content=\"text/html;charset=utf-8\" />";
    ofs << "<title>gencore report at " + getCurrentSystemTime() + " </title>";
    printJS(ofs);
    printCSS(ofs);
    ofs << "</head>";
    ofs << "<body><div id='container'>";
}

void HtmlReporter::printCSS(ofstream& ofs){
    ofs << "<style type=\"text/css\">" << endl;
    ofs << "td {border:1px solid #dddddd;padding:5px;font-size:12px;}" << endl;
    ofs << "table {border:1px solid #999999;padding:2x;border-collapse:collapse; width:800px}" << endl;
    ofs << ".col1 {width:280px; font-weight:bold;}" << endl;
    ofs << ".adapter_col {width:500px; font-size:10px;}" << endl;
    ofs << "img {padding:30px;}" << endl;
    ofs << "#menu {font-family:Consolas, 'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs << "#menu a {color:#0366d6; font-size:18px;font-weight:600;line-height:28px;text-decoration:none;font-family:-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol'}" << endl;
    ofs << "a:visited {color: #999999}" << endl;
    ofs << ".alignleft {text-align:left;}" << endl;
    ofs << ".alignright {text-align:right;}" << endl;
    ofs << ".figure {width:800px;height:600px;}" << endl;
    ofs << ".header {color:#ffffff;padding:1px;height:20px;background:#000000;}" << endl;
    ofs << ".section_title {color:#ffffff;font-size:20px;padding:5px;text-align:left;background:#663355; margin-top:10px;}" << endl;
    ofs << ".subsection_title {font-size:16px;padding:5px;margin-top:10px;text-align:left;color:#663355}" << endl;
    ofs << "#container {text-align:center;padding:3px 3px 3px 10px;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs << ".menu_item {text-align:left;padding-top:5px;font-size:18px;}" << endl;
    ofs << ".highlight {text-align:left;padding-top:30px;padding-bottom:30px;font-size:20px;line-height:35px;}" << endl;
    ofs << "#helper {text-align:left;border:1px dotted #fafafa;color:#777777;font-size:12px;}" << endl;
    ofs << "#footer {text-align:left;padding:15px;color:#ffffff;font-size:10px;background:#663355;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}" << endl;
    ofs << ".kmer_table {text-align:center;font-size:8px;padding:2px;}" << endl;
    ofs << ".kmer_table td{text-align:center;font-size:8px;padding:0px;color:#ffffff}" << endl;
    ofs << ".sub_section_tips {color:#999999;font-size:10px;padding-left:5px;padding-bottom:3px;}" << endl;
    ofs << ".coverage_div {}" << endl;
    ofs << ".bed_coverage_div {}" << endl;
    ofs << "</style>" << endl;
}

void HtmlReporter::printJS(ofstream& ofs){
    ofs << "<script src='https://cdn.plot.ly/plotly-latest.min.js'></script>" << endl;
    ofs << "\n<script type=\"text/javascript\">" << endl;
    ofs << "    function showOrHide(divname) {" << endl;
    ofs << "        div = document.getElementById(divname);" << endl;
    ofs << "        if(div.style.display == 'none')" << endl;
    ofs << "            div.style.display = 'block';" << endl;
    ofs << "        else" << endl;
    ofs << "            div.style.display = 'none';" << endl;
    ofs << "    }" << endl;
    ofs << "</script>" << endl;
}

const string HtmlReporter::getCurrentSystemTime()
{
  auto tt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  struct tm* ptm = localtime(&tt);
  char date[60] = {0};
  sprintf(date, "%d-%02d-%02d      %02d:%02d:%02d",
    (int)ptm->tm_year + 1900,(int)ptm->tm_mon + 1,(int)ptm->tm_mday,
    (int)ptm->tm_hour,(int)ptm->tm_min,(int)ptm->tm_sec);
  return std::string(date);
}

void HtmlReporter::printFooter(ofstream& ofs){
    ofs << "\n</div>" << endl;
    ofs << "<div id='footer'> ";
    ofs << "<p>"<<command<<"</p>";
    ofs << "gencore " << VERSION_NUMBER << ", at " << getCurrentSystemTime() << " </div>";
    ofs << "</body></html>";
}
