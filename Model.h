
#ifndef PinT_MODEL_H
#define PinT_MODEL_H 1

#include "common.h"

/**
 * in current, the model is only for holding the physical variables.
 * In the further, it can be extended to a mesh structure
 */

struct Model {

public:

    //the size of time-space domain
    double tspan = 4.0;
    double xspan = 1.0; 

    void set_tspan(double t){
        this->tspan = t;
    }
    void set_xspan(double x){
        this->xspan = x;
    }

    //the time steps of the whole time domain, fine solver
    long   nt = 200000;

    int   rfc_ = 20;  // fine steps / coarse steps 
    
    void set_nt(long nt){
        this->nt   = nt;
    }
    void set_rfc_(int rfc){
        this->rfc_ = rfc;
    }

    //the number of time slices, parallel processes along time domain 
    int slices = 4;
    void set_slices(int n){
        this->slices = n;
    }

    //the mesh size of the whole space grid
    int  nx = 100;
    void set_nx(int n){
        this->nx = n;
    }

    int nguard = 1;
    void set_nguard(int n){
        this->nguard = n;
    }

    int size;
       
    // the step width of time-space domain
    double dt; 
    double dx;

    int f_steps ;   // number of time steps of fine solver within one time slice 
    int c_steps ;   // number of time steps of coarse solver within one time slice  
    int f_dt  ;     // fine solver
    int c_dt ;      // time step with of coarse solver
  
        
    double *x;

    Model() { }

    void init(){
        size = nx+2*nguard;
        dt = tspan/nt; 
        dx = xspan / (nx+2*nguard-1);

        f_steps = nt/slices ; 
        c_steps = f_steps/rfc_ ;    
        f_dt = dt;              
        c_dt = dt*rfc_; 

        x = alloc_mem(size);

        init_x();
    }

    ~Model(){
        free_mem(x);
    }
   

    virtual int init_x()=0;

    //boundary condition
    virtual void bc()=0;
};

#endif
