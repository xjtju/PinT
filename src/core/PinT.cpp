#include "PinT.h"

PinT* PinT::s_instance = 0;


void PinT::set_tdomain(double tspan, long nt, int rfc, int tsnum) {

    this->Tspan = tspan;
    this->Nt    = nt;
    this->rfc_  = rfc;
    this->tsnum = tsnum;
}

void PinT::set_sdomain(double xspan, long nx, int spnum) {
    this->Xspan = xspan;
    this->Nx    = nx;    
    this->spnum = spnum;
}
    
void PinT::init(){
    f_steps = Nt/tsnum ; 
    c_steps = f_steps/rfc_ ;    

    f_dt = Tspan/Nt;              
    c_dt = f_dt*rfc_; 

    sub_nx = Nx / spnum;
    dx = Xspan / Nx; 

    kpar_limit = tsnum; 
}

void PinT::print() {

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
