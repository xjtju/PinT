#include "PinT.h"

PinT* PinT::s_instance = 0;

// parse ini file and initialize parameters
int PinT::init(char* fname){

    this->ini_file = fname;
    int parse_ret = ini_parse(ini_file, handler, this);  

    if (parse_ret < 0) return -1;

    if(pipelined == 0)
        kpar_limit = tsnum;

    if( test_serial > 0 ) { // for running the serial test, coarse solver->fine solver->coarse solver
        this->tsnum  = 1;
        this->spnumx = 1;
        this->spnumy = 1;
        this->spnumz = 1;
        printf("INFO : SERIAL mode is activated, run parareal iteration only once without any space parallel support.\n"); 
        printf("       and time-space divisions are forcely set to ONE\n\n");
    }

    f_steps = Nt/tsnum ; 
    c_steps = f_steps/rfc_ ;    

    f_dt = Tspan/Nt;              
    c_dt = f_dt*rfc_; 

    nx = Nx / spnumx;
    dx = Xspan / Nx; 

    ny = 1;
    nz = 1;

    if(ndim>=2) { 
        ny = Ny / spnumy;
        dy = Yspan / Ny; 
    }
    if(ndim==3) {
        nz = Nz / spnumz;
        dz = Zspan / Nz; 
    }
    if(ndim==1) spnum = spnumx;
    else if(ndim==2) spnum = spnumx*spnumy;
    else if(ndim==3) spnum = spnumx*spnumy*spnumz;

    return parse_ret;
}

//check the configuration is consist with the real run time
bool PinT::check(){
    bool flag = true;
    
    if(ndim>=4) {
        flag = false;
        if(myid==0) fprintf(stderr, "ERROR : 4D or higher space is not supported.\n");
    }

    if(tsnum != (numprocs/spnum)) {
        flag = false;
        if(myid==0) fprintf(stderr, "ERROR : (the space num)*(the time num) should be equal with the total process num.\n");
    }

    if( (0 != Nx%spnumx) 
        || ( (ndim>=2) && ( 0 != Ny%spnumy) )
        || ( (ndim>=3) && ( 0 != Nz%spnumz) ) ) 
    {
        flag = false; 
        if(myid==0) fprintf(stderr, "WARN  : the cell number is not well divided by the parallel cores.\n");
    }

    if( (test_serial>0) && ( numprocs > 1) ){
        flag = false;
        if(myid==0) fprintf(stderr, "ERROR : SERIAL mode is set, but the total number of CPU cores is %d, greater than ONE.\n", spnum);
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
    printf("  fine   dt        : %e\n", f_dt);
    printf("  coarse dt        : %e\n", c_dt);
    printf("  rfc_ (steps )    : %d\n", rfc_);

    printf("  PIPELINED        : %d\n", pipelined);
    printf("  kpar_limit       : %d\n", kpar_limit);
    printf("  linear solver    : %d\n", linear_solver);
    printf("  converge eps     : %e\n", converge_eps);
    printf("  small residual   : %e\n", smlr);

    printf("  debug out prefix : %s\n", debug_pre);
    printf("  monitor   prefix : %s\n", monitor_pre);
    printf("  coord is output  : %d\n", with_coord);

    printf("  verbose output   : %d\n", verbose);

    printf("  test serial mode : %d\n", test_serial);
    printf("  dump init vars   : %d\n", dump_init);

    printf("\n");

}

int PinT::init_module(void *obj, ini_handler handler) {
    return ini_parse(ini_file, handler, obj);  
}

// the global ini handler
int handler(void* pint, const char* section, const char* name, const char* value)
{
    PinT* conf = (PinT*)pint;

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
    else if (MATCH("parareal", "spnumx")) { conf->spnumx = atoi(value); } 
    else if (MATCH("parareal", "spnumy")) { conf->spnumy = atoi(value); } 
    else if (MATCH("parareal", "spnumz")) { conf->spnumz = atoi(value); } 

    else if (MATCH("parareal", "pipelined"))  { conf->pipelined  = atoi(value); } 
    else if (MATCH("parareal", "kpar_limit")) { conf->kpar_limit = atoi(value); } 
    else if (MATCH("parareal", "rfc_")) { conf->rfc_ = atoi(value); } 

    else if (MATCH("parareal", "linear_solver")){ conf->linear_solver = atoi(value); } 
    else if (MATCH("parareal", "converge_eps")) { conf->converge_eps = atof(value); } 
    else if (MATCH("parareal", "sml_res")) { conf->smlr = atof(value); } 
    
    else if (MATCH("monitor", "debug_pre")) { conf->debug_pre = strdup(value); } 
    else if (MATCH("monitor", "monitor_pre")) { conf->monitor_pre = strdup(value); } 
    else if (MATCH("monitor", "with_coord")) { conf->with_coord = atoi(value); } 
    else if (MATCH("monitor", "verbose")) { conf->verbose = atoi(value); } 
    else if (MATCH("monitor", "test_serial")) { conf->test_serial = atoi(value); } 
    else if (MATCH("monitor", "dump_init")) { conf->dump_init = atoi(value); } 

    else {
        //printf("WARN : unknown ini parameter [%s]/[%s] , ignored. \n", section, name);
        return 1;  
    }
    return 0;
}
