#ifndef PinT_H_
#define PinT_H_ 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "ini.h"

/**
 * the configuration of Parareal method, the current structure is proper when there is only one physical process
 *
 * singleton model
 *
 * thanks to .INI parser https://github.com/benhoyt/inih  
 **/

#define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0

class PinT {

private:
    static PinT* s_instance;
    char* ini_file;

    PinT() {}

public:    

    static PinT* instance()
    {
        if(s_instance == 0)
            s_instance = new PinT();
        return s_instance;
    }

    int init(char* fname);
    void print();
    bool check();      
    // the interface for problem-specific init 
    int init_module(void *obj, ini_handler handler);

    // very important, the settings of x/y/z must be consistent with dims  
    int ndim = 1; 

    //the size of time-space domain
    double Tspan; 
    
    double Xspan; 
    double Yspan = 0; 
    double Zspan = 0; 


    long   Nt;  // the time steps of the whole time domain serially if using fine solver
    long   Nx;  // the whole grid space size covering the whole space domain 
    long   Ny = 1;
    long   Nz = 1;

    int   rfc_ ;  // fine steps / coarse steps 

    int numprocs; // the number of all the CPU cores, dynamically obtained by Driver class.
    int tsnum ; // the number of time partitions(slices), parallel processes along time domain 
    
    int spnum ; // the number of space partitions, parallel processes along the space domain 
    int spnumx;
    int spnumy = 0;
    int spnumz = 0;
    
    /* in a good design , the nx and dx etc. should be in Grid object */
    long nx ; // the subgrid (local) size of one MPI process 
    long ny = 1 ; 
    long nz = 1 ; 
   
    double dx;  // cell width 
    double dy = 1;
    double dz = 1;
    
    int nguard = 1; // the nguard cell number 
    int bc_type =0; // boundary type. 0:always unchanged, default is ZERO; 1: reflect
    double bc_val = 0.0; // valid when bc_type=0 

    int pipelined = 0;
    int kpar_limit;
    int skip_mode = 0;  // default is ZERO, not skip
    // relaxation factor to speed up convergence, 
    // details at (Shulin Wu, 2009, Parareal-Richardson Algorithm) 
    double relax_factor = 1.0 ; 

    // only run the fine solver in serial mode, for measuring the orininal single process/thread performance  
    // NOTE: if the flag is valid, all the parameters related with time-space parallel will be ignored 
    int test_serial = 0; 

    int dump_init = 0; // ouput the initial data before starting   

    // the number of timesteps of fine/coarse solver in one time slice
    long f_steps; 
    long c_steps;
    // time step width of fine/coarse solver 
    double f_dt; 
    double c_dt; 

    // because the following dynamic information will not change during the whole running time
    // for conveniently shareing, the dynamic information is placed in the static configuration class  
    int myid;   // MPI process id in global 
    int mysid;  // process id in space domain
    int mytid;  // process id in time domain (slice number)
    MPI_Comm *sp_comm; //space within the same time slice

    double converge_eps = 1.0e-6;
    double smlr = 1.0e-15;

    int ls_solver = 0;  // default linear solver, 0: BiCG; 1: SOR; ; >1: NULL
    // the parameters of linear solver should be moved to corresponding solver class in later
    // but I wonder it maybe introduce confusion if multiple instances of a solver need different settings.   
    // at current, for debugging and paratmer tunning easily  
    int ls_itmax = 20;      
    double ls_eps = 1.0e-6;
    int ls_abort  = 0;  // when convergence condition is not satisfied until the max iterations, break down the program or not 
    bool ls_precond = false; 
    double ls_relaxfactor = 1.5;

    char* debug_pre;
    char* monitor_pre;
    int with_coord = 1;

    int verbose = 0 ; // output some verbose information for debugging during execution

    
    const static int LS_BiCG_ID= 0;
    const static int LS_SOR_ID = 1;

    const char* LS_SOR = "SOR";
    const char* LS_BiCG= "BiCG";
    
};

int handler(void* pint, const char* section, const char* name, const char* value);

#endif
