#ifndef PinT_DRIVER_H_
#define PinT_DRIVER_H_ 1

#include <stdio.h>
#include <stdarg.h>
#include <mpi.h>

#include "ini.h"
#include "PinT.h"
#include "blas.h"
#include "Solver.h"
#include "Grid.h"

/**
 * The driver of the PinT process
 * it implements the Parareal algorithm, and controls the execution flow of the program, especially the time parallel flow.
 * 
 * the core function of the Driver is evolve(), it provides the template of Parareal algorithm, 
 * and drives the problem-specific coarse/fine solvers to evolve along the time slices within the whole time domain.    
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

    void evolve(Grid* g, Solver* G, Solver* F); // the algorithm framework of Parareal method

    void finalize();

    void INFO (const char* fmt, ...);
    void WARN (const char* fmt, ...);
    void Abort(const char* fmt, ...);
  
    double smlr = 1.0e-12;

    /**
     * Parareal's iterative formula : F = G + F - G  
     *
     * NOTE:
     * Due to the "Machine Epsilon" or rounding error, residual calculation is very important.
     * Sometimes, in theory, the residual should be ZERO, but in practice, the calculation value is not ZERO, 
     * despite it is very small, it will has an unignorable impact on convergency due to  accumulating effect.
     * the "smlr" is the threshold for residual control, when res < smlr, res is regarded as ZERO.
     */
    inline void pint_sum(Grid *grid, double *u, double *f, double *g, double *g_, double *res, double *sml)  {
        switch(grid->ndim) {
            case 1: blas_pint_sum_1d_(grid->nxyz, &grid->nguard, u, f, g, g_, res, sml); break;
            case 2: blas_pint_sum_2d_(grid->nxyz, &grid->nguard, u, f, g, g_, res, sml); break;
            case 3: printf("ERROR: 3D is not finished"); break;
        }
    }
    
    // check the residual is whether enough small as the Kpar increasing, for DEBUG only 
    void monitorResidual(Grid *g, double res_loc, double max_res,int size );
    // used for residual check
    inline void vector_dist(Grid *g, double *d, double *s, double *val) {
        switch(g->ndim) {
            case 1: blas_vdist_1d_(g->nxyz, &g->nguard, d, s, val); break;
            case 2: blas_vdist_2d_(g->nxyz, &g->nguard, d, s, val); break;
            case 3: printf("ERROR: 3D is not finished"); break;
        }
    }

    /**
     * check whether the current slice is the first, 
     * if it is the first time slice, it should not receive result from  any other slices.     
     * when there is only one MPI process, it is the first slice and also the last slice
     */
    inline bool isFirstSlice(int myid) {
        return (mytid == 0) ? true : false ; 
    }

    /**
     * check whether the current slice is the last, 
     * if it is the last time slice, it should not send result to any other slices.     
     */
    inline bool isLastSlice(int myid) {
        return (mytid == (tsnum-1)) ? true : false ;
    }

    inline int getSliceNum(int myid) {
        return mytid;
    }

};

#endif
