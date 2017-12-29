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
    static PinT* s_instance ;

    PinT() {}

public:    

    static PinT* instance()
    {
        if(s_instance == 0 )
            s_instance = new PinT();
        return s_instance;
    }

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
    
   
    void set_tdomain(double tspan, long nt, int rfc, int tsnum) {

        this->Tspan = tspan;
        this->Nt    = nt;
        this->rfc_  = rfc;
        this->tsnum = tsnum;
    }

    void set_sdomain(double xspan, long nx, int spnum) {
        this->Xspan = xspan;
        this->Nx    = nx;    
        this->spnum = spnum;
    }
    
    void init(){
        f_steps = Nt/tsnum ; 
        c_steps = f_steps/rfc_ ;    

        f_dt = Tspan/Nt;              
        c_dt = f_dt*rfc_; 

        sub_nx = Nx / spnum;
        dx = Xspan / Nx; 

        kpar_limit = tsnum; 
    }

    void print() {

        printf("PinT ini configuration : \n");
        printf("  space dimension  : %d\n", dims);
        printf("  space domain     : [%f, %f, %f]\n", Xspan, Yspan, Zspan);
        printf("  grid cell        : [%f, %f, %f]\n", dx, dx,dx);
        printf("  sub space domain : [%d, %d, %d]\n", sub_nx, sub_nx, sub_nx);
        printf("  guard cells      : %d\n",nguard);

        printf("  time domin       : %f\n", Tspan);

        printf("  space parallel   : %d\n", spnum);
        printf("  time  parallel   : %d\n", tsnum);

        printf("  serial time steps: %d\n", Nt);
        printf("  fine   steps     : %d\n", f_steps);
        printf("  coarse steps     : %d\n", c_steps);
        printf("  fine   dt        : %f\n", f_dt);
        printf("  coarse dt        : %f\n", c_dt);
        printf("  rfc_ (steps )    : %d\n", rfc_);
        printf("  kpar_limit       : %d\n", kpar_limit);
        printf("  converge eps     : %f\n", converge_eps);
        printf("\n");
    }
};
PinT* PinT::s_instance = 0;
#endif


