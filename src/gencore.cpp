#include "gencore.h"
#include "bamutil.h"
#include "jsonreporter.h"
#include "htmlreporter.h"

Gencore::Gencore(Options *opt){
    mOptions = opt;
    mBamHeader = NULL;
    mOutSam = NULL;
    mPreStats = new Stats(opt);
    mPostStats = new Stats(opt);
}

Gencore::~Gencore(){
    releaseClusters(mProperClusters);
    releaseClusters(mUnProperClusters);
    if(mBamHeader != NULL) {
        bam_hdr_destroy(mBamHeader);
        mBamHeader = NULL;
    }
    if(mOutSam != NULL) {
        if (sam_close(mOutSam) < 0) {
            cerr << "ERROR: failed to close " << mOutput << endl;
            exit(-1);
        }
    }
    delete mPreStats;
    delete mPostStats;
}

void Gencore::report() {
    JsonReporter jsonreporter(mOptions);
    jsonreporter.report(mPreStats, mPostStats);
    HtmlReporter htmlreporter(mOptions);
    htmlreporter.report(mPreStats, mPostStats);
}

void Gencore::releaseClusters(map<int, map<int, map<int, Cluster*>>>& clusters) {
    map<int, map<int, map<int, Cluster*>>>::iterator iter1;
    map<int, map<int, Cluster*>>::iterator iter2;
    map<int, Cluster*>::iterator iter3;
    for(iter1 = clusters.begin(); iter1 != clusters.end(); iter1++) {
        for(iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++) {
            for(iter3 = iter2->second.begin(); iter3 != iter2->second.end(); iter3++) {
                delete iter3->second;
            }
        }
    }
}

void Gencore::dumpClusters(map<int, map<int, map<int, Cluster*>>>& clusters) {
    map<int, map<int, map<int, Cluster*>>>::iterator iter1;
    map<int, map<int, Cluster*>>::iterator iter2;
    map<int, Cluster*>::iterator iter3;
    for(iter1 = clusters.begin(); iter1 != clusters.end(); iter1++) {
        for(iter2 = iter1->second.begin(); iter2 != iter1->second.end(); iter2++) {
            for(iter3 = iter2->second.begin(); iter3 != iter2->second.end(); iter3++) {
                 iter3->second->dump();
            }
        }
    }
}

void Gencore::outputPair(Pair* p) {
    if(mOutSam == NULL || mBamHeader == NULL)
        return ;

    if(p->mLeft) {
        mPostStats->addRead(p->mLeft->core.l_qseq, BamUtil::getED(p->mLeft));
        mPostStats->statDepth(p->mLeft->core.tid, p->mLeft->core.pos, p->mLeft->core.l_qseq);
        if(sam_write1(mOutSam, mBamHeader, p->mLeft) <0) {
            cerr << "Writing failed, exiting..." << endl;
            exit(-1);
        }
    }
    if(p->mRight) {
        mPostStats->addRead(p->mRight->core.l_qseq, BamUtil::getED(p->mRight));
        mPostStats->statDepth(p->mRight->core.tid, p->mRight->core.pos, p->mRight->core.l_qseq);
        if(sam_write1(mOutSam, mBamHeader, p->mRight) <0) {
            cerr << "Writing failed, exiting..." << endl;
            exit(-1);
        }
    }
    mPostStats->addMolecule(1, p->mLeft && p->mRight);
}

void Gencore::consensus(){
    samFile *in;
    in = sam_open(mOptions->input.c_str(), "r");
    if (!in) {
        cerr << "ERROR: failed to open " << mOptions->input << endl;
        exit(-1);
    }

    if(ends_with(mOptions->output, "sam"))
        mOutSam = sam_open(mOptions->output.c_str(), "w");
    else 
        mOutSam = sam_open(mOptions->output.c_str(), "wb");
    if (!mOutSam) {
        cerr << "ERROR: failed to open output " << mOptions->output << endl;
        exit(-1);
    }

    mBamHeader = sam_hdr_read(in);
    mOptions->bamHeader = mBamHeader;
    mPreStats->makeGenomeDepthBuf();
    mPostStats->makeGenomeDepthBuf();
    if (mBamHeader == NULL || mBamHeader->n_targets == 0) {
        cerr << "ERROR: this SAM file has no header " << mInput << endl;
        exit(-1);
    }
    BamUtil::dumpHeader(mBamHeader);

    if (sam_hdr_write(mOutSam, mBamHeader) < 0) {
        cerr << "failed to write header" << endl;
        exit(-1);
    }

    bam1_t *b = NULL;
    b = bam_init1();
    int r;
    int count = 0;
    int lastTid = -1;
    int lastPos = -1;
    while ((r = sam_read1(in, mBamHeader, b)) >= 0) {
        // check whether the BAM is sorted
        if(b->core.tid <lastTid || (b->core.tid == lastTid && b->core.pos <lastPos)) {
            // skip the -1:-1, which means unmapped
            if(b->core.tid >=0 && b->core.pos >= 0) {
                cerr << "ERROR: the input is unsorted. Found unsorted read in " << b->core.tid << ":" << b->core.pos << endl;
                cerr << "Please sort the input first." << endl << endl;
                exit(-1);
            }
        }
        // for testing, we only process to some contig
        if(mOptions->maxContig>0 && b->core.tid>=mOptions->maxContig){
            b = bam_init1();
            break;
        }
        // if debug flag is enabled, show which contig we are start to process
        if(mOptions->debug && b->core.tid > lastTid) {
            cerr << "Starting contig " << b->core.tid << endl;
        }
        lastTid = b->core.tid;
        lastPos = b->core.pos;

        // unmapped reads, we just write it and continue
        if(b->core.tid < 0 || b->core.pos < 0) {
            mPreStats->addRead(b->core.l_qseq, 0, false);
            mPostStats->addRead(b->core.l_qseq, 0, false);
            if(sam_write1(mOutSam, mBamHeader, b) <0) {
                cerr << "Writing failed, exiting..." << endl;
                exit(-1);
            }
            continue;
        }

        // for secondary alignments, we just skip it
        if(!BamUtil::isPrimary(b)) {
            continue;
        }
        mPreStats->statDepth(b->core.tid, b->core.pos, b->core.l_qseq);
        addToCluster(b);
        b = bam_init1();
    }

    finishConsensus(mProperClusters);
    finishConsensus(mUnProperClusters);

    bam_destroy1(b);
    sam_close(in);

    cerr << "----Before gencore processing:" << endl;
    mPreStats->print();

    cerr << endl << "----After gencore processing:" << endl;
    mPostStats->print();

    report();
}

void Gencore::addToProperCluster(bam1_t* b) {
    int tid = b->core.tid;
    int tlen = b->core.isize;

    // add a WAR fix for some unproperly mapped pairs
    /*if(abs(b->core.isize) > 10000) {
        tlen = b->core.mpos - b->core.pos;
    }*/

    int left = b->core.pos;
    if(tlen < 0)
        left = b->core.mpos;
    int right = left + abs(tlen) - 1;
    createCluster(mProperClusters, tid, left, right);
    mProperClusters[tid][left][right]->addRead(b);


    static int tick = 0;
    tick++;
    if(tick % 10000 != 0)
        return;

    // make consensus merge
    map<int, map<int, map<int, Cluster*>>>::iterator iter1;
    map<int, map<int, Cluster*>>::iterator iter2;
    map<int, Cluster*>::iterator iter3;
    bool needBreak = false;
    for(iter1 = mProperClusters.begin(); iter1 != mProperClusters.end();) {
        if(iter1->first > tid || needBreak)
            break;
        for(iter2 = iter1->second.begin(); iter2 != iter1->second.end(); ) {
            if(iter1->first == tid && iter2->first >= b->core.pos) {
                needBreak = true;
                break;
            }
            for(iter3 = iter2->second.begin(); iter3 != iter2->second.end(); ) {
                // only deal with the clusters with right < processing pos
                if(iter1->first == tid && iter3->first >= b->core.pos) {
                    break;
                }
                vector<Pair*> csPairs = iter3->second->clusterByUMI(mOptions->properReadsUmiDiffThreshold, mPreStats, mPostStats);
                for(int i=0; i<csPairs.size(); i++) {
                    //csPairs[i]->dump();
                    outputPair(csPairs[i]);
                    delete csPairs[i];
                }
                // this tid:left:right is done
                delete iter3->second;
                iter3 = iter2->second.erase(iter3);
            }
            // this tid:left is done
            if(iter2->second.size() == 0) {
                iter2 = iter1->second.erase(iter2);
            } else {
                iter2++;
            }
        }
        // this tid is done
        if(iter1->second.size() == 0) {
            iter1 = mProperClusters.erase(iter1);
        } else {
            iter1++;
        }
    }
}

void Gencore::finishConsensus(map<int, map<int, map<int, Cluster*>>>& clusters) {
    // make consensus merge
    map<int, map<int, map<int, Cluster*>>>::iterator iter1;
    map<int, map<int, Cluster*>>::iterator iter2;
    map<int, Cluster*>::iterator iter3;
    for(iter1 = clusters.begin(); iter1 != clusters.end();) {
        for(iter2 = iter1->second.begin(); iter2 != iter1->second.end(); ) {
            for(iter3 = iter2->second.begin(); iter3 != iter2->second.end(); ) {
                // for unmapped reads, we just store them
                if(iter1->first < 0 || iter2->first < 0 || iter3->first < 0 ) {
                    map<string, Pair*>::iterator iterOfPairs;
                    for(iterOfPairs = iter3->second->mPairs.begin(); iterOfPairs!=iter3->second->mPairs.end(); iterOfPairs++) {
                        //csPairs[i]->dump();
                        outputPair(iterOfPairs->second);
                    }
                } else {
                    vector<Pair*> csPairs = iter3->second->clusterByUMI(mOptions->unproperReadsUmiDiffThreshold, mPreStats, mPostStats);
                    for(int i=0; i<csPairs.size(); i++) {
                        //csPairs[i]->dump();
                        outputPair(csPairs[i]);
                        delete csPairs[i];
                    }
                }
                // this tid:left:right is done
                delete iter3->second;
                iter3 = iter2->second.erase(iter3);
            }
            // this tid:left is done
            if(iter2->second.size() == 0) {
                iter2 = iter1->second.erase(iter2);
            } else {
                iter2++;
            }
        }
        // this tid is done
        if(iter1->second.size() == 0) {
            iter1 = clusters.erase(iter1);
        } else {
            iter1++;
        }
    }
}

void Gencore::addToUnProperCluster(bam1_t* b) {
    int tid = b->core.tid;
    int left = b->core.pos;
    int right = b->core.mpos;
    if(b->core.mtid < b->core.tid) {
        tid = b->core.mtid;
        left = b->core.mpos;
        right = b->core.pos;
    }
    createCluster(mUnProperClusters, tid, left, right);
    mUnProperClusters[tid][left][right]->addRead(b);
}

void Gencore::createCluster(map<int, map<int, map<int, Cluster*>>>& clusters, int tid, int left, int right) {
    map<int, map<int, map<int, Cluster*>>>::iterator iter1 = clusters.find(tid);

    if(iter1 == clusters.end()) {
        clusters[tid] = map<int, map<int, Cluster*>>();
        clusters[tid][left] = map<int, Cluster*>();
        clusters[tid][left][right] = new Cluster(mOptions);
    } else {
        map<int, map<int, Cluster*>>::iterator iter2  =iter1->second.find(left);
        if(iter2 == iter1->second.end()) {
            clusters[tid][left] = map<int, Cluster*>();
            clusters[tid][left][right] = new Cluster(mOptions);
        } else {
            map<int, Cluster*>::iterator iter3 = iter2->second.find(right);
            if(iter3 == iter2->second.end())
                clusters[tid][left][right] = new Cluster(mOptions);
        }
    }
}

void Gencore::addToCluster(bam1_t* b) {
    mPreStats->addRead(b->core.l_qseq, BamUtil::getED(b));
    // unproperly mapped
    if(b->core.tid < 0) {
        addToUnProperCluster(b);
    } else {
        addToProperCluster(b);
    }
}