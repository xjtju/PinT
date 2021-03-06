#include <stdio.h>
#include <unistd.h> 
#include "Monitor.h"

int Monitor::counter = -1;

std::string Monitor::CSolver = "C Solver";
std::string Monitor::FSolver = "F Solver";
std::string Monitor::SEND    = "Send Slice";
std::string Monitor::RECV    = "Recv Slice";
std::string Monitor::GC      = "Guard Cell";
std::string Monitor::RES     = "Residual";

#ifdef _PMLib_
PerfMonitor Monitor::pm = PerfMonitor();
#endif

void Monitor::initLabels() {
    if(counter>=0) return;

    setLabel_CALC(CSolver, true);
    setLabel_CALC(FSolver, true);
    setLabel_CALC(RES, true);

    setLabel_COMM(SEND, true);
    setLabel_COMM(RECV, true);

    setLabel_COMM(GC,   true);
}

// the wraper functions of PMLib
void Monitor::initialize() {
#ifdef _PMLib_
    pm.initialize();  
    initLabels();
    counter = 0;
#endif
}

void Monitor::setRankInfo(int id) {
#ifdef _PMLib_
    pm.setRankInfo(id);
    int num_thread = omp_get_max_threads();
    if(id==0) printf("INFO : PMLib is used for performance monitor: %d threads per process!\n", num_thread);
#else
    if(id==0) printf("INFO : Performance monitor is not configurated, no performance data will be collected !\n");
#endif 
}

void Monitor::setLabel_CALC(const std::string &label, bool exclusive) {
#ifdef _PMLib_
    pm.setProperties(label, PerfMonitor::CALC, exclusive);
    counter ++;
#endif 
}

void Monitor::setLabel_COMM(const std::string &label, bool exclusive) {
#ifdef _PMLib_
    pm.setProperties(label, PerfMonitor::COMM, exclusive);
    counter ++;
#endif 
}

void Monitor::start(const std::string &label){
#ifdef _PMLib_
    pm.start(label);
#endif 
}
void Monitor::stop (const std::string &label, double flopPerTask, unsigned iterationCount){
#ifdef _PMLib_
    pm.stop(label,flopPerTask,iterationCount);
#endif 
}

void Monitor::gather(){
#ifdef _PMLib_
    pm.gather();
#endif 
}

void Monitor::print(FILE* fp, const std::string hostname, const std::string comments, int seqSections) {
#ifdef _PMLib_
    pm.print(fp, hostname, comments, seqSections);
#endif 
}
void Monitor::printDetail(FILE* fp, int legend, int seqSections) {
#ifdef _PMLib_
    pm.printDetail(fp, legend, seqSections);
#endif 
}


/**
 * Form "man 2 close" : 
 *  A successful close does not guarantee that the data has been successfully saved to disk, as the kernel defers writes.
 *  The man page says that if you want to be sure that your data are on disk, you have to use fsync() yourself.
 *
 * THOUGH I had used all the tips, it still cannot make sure that the profiling information is completely written to the disk on the ITO supercomputer of Kyusyu university. But if using stdout/stderr instead of common files, it is OK. 
 */
void Monitor::print(const char* fname, const std::string hostname, const std::string comments, int seqSections){
#ifdef _PMLib_
    FILE *fp;
    fp = fopen(fname,"w");
    pm.print(fp, hostname, comments, seqSections);
    fflush(fp);
   // int fd = fileno(fp);
   // fsync(fd);
    fclose(fp);
#endif
}

void Monitor::printDetail(const char* fname, int legend, int seqSections){
#ifdef _PMLib_
    FILE *fp;
    fp = fopen(fname,"w");
    pm.printDetail(fp, legend, seqSections);
    fflush(fp);
    //int fd = fileno(fp);
    //fsync(fd);
    fclose(fp);
#endif
}
