#include "unittest.h"
#include "bamutil.h"
#include <time.h>
#include "cluster.h"

UnitTest::UnitTest(){

}

void UnitTest::run(){
    bool passed = true;
    passed &= BamUtil::test();
    passed &= Cluster::test();
    printf("\n==========================\n");
    printf("%s\n\n", passed?"PASSED":"FAILED");
}