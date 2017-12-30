#ifndef PinT_H_
#define PinT_H_ 1

#include <mpi.h>
#include <stdio.h>

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
    void set_tdomain(double tspan, long nt, int rfc, int tsnum) ;
    void set_sdomain(double xspan, long nx, int spnum) ;

public:    

    static PinT* instance()
    {
        if(s_instance == 0 )
            s_instance = new PinT();
        return s_instance;
    }

    void init();
    void print() ;

    //the size of time-space domain
    double Tspan; 

    int dims;
    double Xspan; 
    double Yspan; 
    double Zspan; 


    long   Nt;  // the time steps of the whole time domain serially if using fine solver
    long   Nx;  // the whole grid space size 
    long   Ny;
    long   Nz;


    int   rfc_ ;  // fine steps / coarse steps 

    int tsnum  ; // the number of time partitions(slices), parallel processes along time domain 
    int spnum  ; // the number of space partitions, parallel processes along the space domain 

    int kpar_limit;

    // the number of timesteps of fine/coarse solver in one time slice
    long f_steps; 
    long c_steps;
    // time step width of fine/coarse solver 
    double f_dt; 
    double c_dt; 
    double dx;  // cell width 

    int sub_nx; // the subgrid size of one MPI process 
    int nguard = 1; // the nguard cell number 

    double converge_eps = 1.0e-6;


};

#endif
