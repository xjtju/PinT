#ifndef PinT_DRIVER_H_
#define PinT_DRIVER_H_ 1

#include <stdio.h>
#include <stdarg.h>
#include <mpi.h>

#include "ini.h"
#include "PinT.h"
#include "Solver.h"
#include "Grid.h"

/**
 * the driver of the PinT process 
 */

class Driver {

    PinT *conf;
    MPI_Comm sp_comm; //space within the same time slice

 public:
    
    int myid = 0;
    int numprocs; // MPI process number
    int spnum ; // space parallel number
    int tsnum ; // time  parallel number

    int mytid ; // time slice number 
    int mysid ; // space subgrid number
    
    int kpar = 0;

    void init(int argc, char* argv[]);

    void evolve(Grid* g, Solver* G, Solver* F);

    void finalize();

    void INFO (const char* fmt, ...);
    void WARN (const char* fmt, ...);
    void Abort(const char* fmt, ...);

    void monitorResidual(double* u_c, double* u_prev, double* u_end, double* u_f, double res_loc, double max_res,int size );

    /**
     * check whether the current slice is the first, 
     * if it is the first time slice, it should not receive result from  any other slices.     
     */
    inline bool isFirstSlice(int myid) {
        return (mytid == 0) ? true : false ; 
    }
    /**
     * check whether the current slice is the last, 
     * if it is the last time slice, it should not send result to any other slices.     
     */
    inline bool isLastSlice(int myid) {
        return (myid/spnum == (tsnum-1)) ? true : false ;
    }

    inline int getSliceNum(int myid) {
        return mytid;
    }

};

#endif
