#include "PinT.h"

PinT* PinT::s_instance = 0;

void PinT::init(){
    f_steps = Nt/tsnum ; 
    c_steps = f_steps/rfc_ ;    

    f_dt = Tspan/Nt;              
    c_dt = f_dt*rfc_; 

    nx = Nx / spnumx;
    dx = Xspan / Nx; 
    

    if(ndim>=2) { 
        ny = Ny / spnumy;
        dy = Yspan / Ny; 
    }
    if(ndim==3) {
        nz = Nz / spnumz;
        dz = Zspan / Nz; 
    }

    kpar_limit = tsnum;
}

//check the configuration is consist with the real run time
bool PinT::check(){
    bool flag = true;
    
    if(ndim>=4) {
        flag = false;
        fprintf(stderr, "ERROR : 4D or higher space is not supported.\n");
    }

    if(tsnum != (numprocs/spnum)) {
        flag = false;
        fprintf(stderr, "ERROR : (the space num)*(the time num) should be equal with the total process num.\n");
    }

    if( spnum != (spnumx*spnumy*spnumz) ) {
        flag = false;
        fprintf(stderr, "ERROR : (the space num) should be equal with value multiplied from all directions.\n");
    }

    if( (0 != Nx%spnumx) 
        || ( (ndim>=2) && ( 0 != Ny%spnumy) )
        || ( (ndim>=3) && ( 0 != Nz%spnumz) ) ) 
    {
        flag = false; 
        fprintf(stderr, "WARN  : the cell number is not well divided by the parallel cores.\n");
    }

    return flag;
}

void PinT::print() {

    printf("PinT ini configuration : \n");
    printf("  space dimensions : %d\n", ndim);

    printf("  space domain size: [%f, %f, %f]\n", Xspan, Yspan, Zspan);
    printf("  unit cell size   : [%f, %f, %f]\n", dx, dy,dz);
    printf("  all space cells  : [%d, %d, %d]\n", Nx, Ny, Nz);
    printf("  sub space cells  : [%d, %d, %d]\n", nx, ny, nz);
    printf("  space division   : [%d, %d, %d]\n", spnumx, spnumy, spnumz);
    printf("  guard cells      : %d\n",nguard);

    printf("  boundary type    : %d\n",bc_type);
    printf("  bc_val(type=0)   : %f\n",bc_val);

    printf("  time domin       : %f\n", Tspan);
    printf("  time  parallel   : %d\n", tsnum);
    printf("  space parallel   : %d\n", spnum);

    printf("  serial time steps: %d\n", Nt);
    printf("  fine   steps     : %d\n", f_steps);
    printf("  coarse steps     : %d\n", c_steps);
    printf("  fine   dt        : %f\n", f_dt);
    printf("  coarse dt        : %f\n", c_dt);
    printf("  rfc_ (steps )    : %d\n", rfc_);
    printf("  kpar_limit       : %d\n", kpar_limit);
    printf("  converge eps     : %e\n", converge_eps);
    printf("  small residual   : %e\n", smlr);

    printf("  debug out prefix : %s\n", debug_pre);
    printf("  monitor   prefix : %s\n", monitor_pre);

    printf("\n");
}

int handler(void* pint, const char* section, const char* name, const char* value)
{
    PinT* conf = (PinT*)pint;

    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
    if (MATCH("domain","ndim"))        { conf->ndim = atoi(value); }

    else if (MATCH("domain", "Tspan")) { conf->Tspan = atof(value); } 
    else if (MATCH("domain", "Nt"))    { conf->Nt = atoi(value); } 

    else if (MATCH("domain", "Xspan")) { conf->Xspan = atof(value); } 
    else if (MATCH("domain", "Yspan")) { conf->Yspan = atof(value); } 
    else if (MATCH("domain", "Zspan")) { conf->Zspan = atof(value); } 

    else if (MATCH("domain", "Nx"))    { conf->Nx = atoi(value); } 
    else if (MATCH("domain", "Ny"))    { conf->Ny = atoi(value); } 
    else if (MATCH("domain", "Nz"))    { conf->Nz = atoi(value); } 
    else if (MATCH("domain", "nguard")){ conf->nguard = atoi(value); } 

    else if (MATCH("domain", "bc_type")){ conf->bc_type = atoi(value); } 
    else if (MATCH("domain", "bc_val")) { conf->bc_val  = atof(value); } 

    else if (MATCH("parareal", "tsnum")) { conf->tsnum = atoi(value); } 
    else if (MATCH("parareal", "spnum")) { conf->spnum = atoi(value); } 
    else if (MATCH("parareal", "spnumx")) { conf->spnumx = atoi(value); } 
    else if (MATCH("parareal", "spnumy")) { conf->spnumy = atoi(value); } 
    else if (MATCH("parareal", "spnumz")) { conf->spnumz = atoi(value); } 

    else if (MATCH("parareal", "kpar_limit")) { conf->kpar_limit = atoi(value); } 
    else if (MATCH("parareal", "rfc_")) { conf->rfc_ = atoi(value); } 

    else if (MATCH("parareal", "converge_eps")) { conf->converge_eps = atof(value); } 
    else if (MATCH("parareal", "sml_res")) { conf->smlr = atof(value); } 
    
    else if (MATCH("monitor", "debug_pre")) { conf->debug_pre = strdup(value); } 
    else if (MATCH("monitor", "monitor_pre")) { conf->monitor_pre = strdup(value); } 

    else {
        printf("WARN : unknown ini parameter [%s]/[%s] , ignored. \n", section, name);
        return 1;  
    }
    return 0;
}
