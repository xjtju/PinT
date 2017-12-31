#ifndef PinT_H_
#define PinT_H_ 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * the configuration of Parareal method, the current structure is proper when there is only one physical process
 *
 * singleton model
 *
 **/
class PinT {

private:
    static PinT* s_instance;

    PinT() {}

public:    

    static PinT* instance()
    {
        if(s_instance == 0 )
            s_instance = new PinT();
        return s_instance;
    }

    void init();
    void print();
      
    // very important, the settings of x/y/z must be consistent with dims  
    int dims = 1; 

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

    int tsnum ; // the number of time partitions(slices), parallel processes along time domain 

    int spnum ; // the number of space partitions, parallel processes along the space domain 
    int spnumx;
    int spnumy = 0;
    int spnumz = 0;

    long nx ; // the subgrid (local) size of one MPI process 
    long ny = 1 ; 
    long nz = 1 ; 

    double dx;  // cell width 
    double dy = 1;
    double dz = 1;

    int nguard = 1; // the nguard cell number 

    int kpar_limit;

    // the number of timesteps of fine/coarse solver in one time slice
    long f_steps; 
    long c_steps;
    // time step width of fine/coarse solver 
    double f_dt; 
    double c_dt; 

    double converge_eps = 1.0e-6;

};

int handler(void* pint, const char* section, const char* name, const char* value);

#endif
