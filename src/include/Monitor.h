#ifndef PinT_Monitor_H_
#define PinT_Monitor_H_ 1


#ifdef _PMLib_
#include <PerfMonitor.h>
using namespace pm_lib;
#endif

#include <string>
using namespace std;

/**
 * the simple wraper of PMLib
 */
class Monitor {

protected:
    static int counter; 

#ifdef _PMLib_
    static PerfMonitor pm;
#endif

public:

   void initialize();
   void setRankInfo(int id);

   void setLabel_CALC(const std::string &label, bool exclusive=true);
   void setLabel_COMM(const std::string &label, bool exclusive=true);

   void start(const std::string &label);
   void stop (const std::string &label, double flopPerTask=0.0, unsigned iterationCount=1);

   void gather();

   void print(FILE* fp, const std::string hostname, const std::string comments, int seqSections=0);
   void printDetail(FILE* fp, int legend=0, int seqSections=0);

   void print(const char* fname, const std::string hostname, const std::string comments, int seqSections=0);
   void printDetail(const char* fname, int legend=0, int seqSections=0);

   // set all the labels monitored in the program 
   void initLabels();

   static std::string RECV ;
   static std::string SEND ;
   static std::string CSolver;
   static std::string FSolver;
   static std::string GC;    // guard cell
   static std::string RES;   //residual
};

#endif
