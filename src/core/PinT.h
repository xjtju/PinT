#ifndef PinT_H_
#define PinT_H_ 1

#include <stdio.h>
/**
 * the configuration of Parareal method, the current structure is proper when there is only one physical process
 *
 * singleton model
 *
 **/
struct PinT {

private:
    PinT() {};
    //PinT(PinT const&);              // Don't Implement
    void operator=(PinT const&);    // Don't implement

public:    

    static PinT& instance()
    {
        static PinT instance; 
        return instance;
    }

    //the size of time-space domain
    double Tspan = 1.5;
    double Xspan = 1.0; 

    inline void set_Tspan(double t){ this->Tspan = t; }
    inline void set_Xspan(double x){ this->Xspan = x; }

    long   Nt = 1000000;  // the time steps of the whole time domain serially if using fine solver
    long   Nx = 100;     // the whole grid space size 
    inline void set_Nt(long nt){ this->Nt = nt; }
    inline void set_Nx(long nx){ this->Nx = nx; }

    int   rfc_ = 50;  // fine steps / coarse steps 
    inline void set_rfc_(int rfc){ this->rfc_ = rfc; }

    int slices    = 4; // the number of time partitions(slices), parallel processes along time domain 
    int subgrids  = 1; // the number of space partitions, parallel processes along the space domain 
    inline void set_slices(int n){ this->slices = n; }
    inline void set_subgrids(int n){ this->subgrids = n; }

    int kpar_limit;
    inline void set_kpar_limit(int k) { this->kpar_limit = k;}
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
    
    
    void set_tdomain(double tspan, long nt, int rfc, int slices) {

        this->Tspan = tspan;
        this->Nt    = nt;
        this->rfc_  = rfc;
        this->slices = slices;
    }

    void set_sdomain(double xspan, long nx, int sgrids) {
        this->Xspan = xspan;
        this->Nx    = nx;    
        this->subgrids = sgrids;
    }
    
    void init(){
        f_steps = Nt/slices ; 
        c_steps = f_steps/rfc_ ;    

        f_dt = Tspan/Nt;              
        c_dt = f_dt*rfc_; 

        sub_nx = Nx / subgrids;
        dx = Xspan / Nx; 

        kpar_limit = slices; 
    }

    void print() {
        printf("PinT configuration : \n");
        printf("  space num : %d\n", subgrids);
        printf("  time  num : %d\n", slices);
        printf("  fine steps: %d\n", f_steps);
        printf("  coar steps: %d\n", c_steps);
        printf("  dx        : %f\n", dx);
        printf("  fine   dt : %f\n", f_dt);
        printf("  coarse dt : %f\n", c_dt);
        printf("  rfc_      : %d\n", rfc_);
        printf("\n");
    }
};
#endif
