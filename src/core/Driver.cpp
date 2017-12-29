#include "Driver.h"


void Driver::init(char* ini_file, PinT* conf){

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    INFO("Initialization file is  %s \n", ini_file);

    if (ini_parse(ini_file, handler, conf) < 0) {
        Abort("Can't load ini file : %s .\n", ini_file);
    } else {
        conf->init();
        if(myid==0) conf->print();
    }
}


void Driver::evolve(){}

void Driver::finalize() {}

void Driver::INFO(const char* fmt, ...) {

    if(myid != 0) return ;
     
    int ret;
    va_list args;

    va_start(args, fmt);
    fputs("INFO : ", stdout);
    ret = vfprintf(stdout, fmt, args);
    va_end(args);
}

void Driver::WARN(const char* fmt, ...) {

    if(myid != 0) return ;
     
    int ret;
    va_list args;

    va_start(args, fmt);
    fputs("WARN : ", stdout);
    ret = vfprintf(stdout, fmt, args);
    va_end(args);
}


void Driver::Abort(const char* fmt, ...) {

    va_list args;
    
    va_start(args, fmt);
    fputs("ERROR: ", stderr);
    vfprintf(stderr, fmt, args);
    va_end(args);

    MPI_Finalize();

    exit(1);
}

static int handler(void* pint, const char* section, const char* name, const char* value)
{
    PinT* conf = (PinT*)pint;

    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
    if (MATCH("domain","dims"))        { conf->dims = atoi(value); }

    else if (MATCH("domain", "Tspan")) { conf->Tspan = atof(value); } 
    else if (MATCH("domain", "Nt"))    { conf->Nt = atoi(value); } 

    else if (MATCH("domain", "Xspan")) { conf->Xspan = atof(value); } 
    else if (MATCH("domain", "Yspan")) { conf->Yspan = atof(value); } 
    else if (MATCH("domain", "Zspan")) { conf->Zspan = atof(value); } 

    else if (MATCH("domain", "Nx"))    { conf->Nx = atoi(value); } 
    else if (MATCH("domain", "Ny"))    { conf->Ny = atoi(value); } 
    else if (MATCH("domain", "Nz"))    { conf->Nz = atoi(value); } 
    else if (MATCH("domain", "nguard")){ conf->nguard = atoi(value); } 

    else if (MATCH("parareal", "tsnum")) { conf->tsnum = atoi(value); } 
    else if (MATCH("parareal", "spnum")) { conf->spnum = atoi(value); } 

    else if (MATCH("parareal", "kpar_limit")) { conf->kpar_limit = atoi(value); } 
    else if (MATCH("parareal", "rfc_")) { conf->rfc_ = atoi(value); } 

    else if (MATCH("parareal", "converge_eps")) { conf->converge_eps = atof(value); } 

    else {
        printf("WARN : unknown ini parameter [%s]/[%s] , ignored. \n", section, name);
        return 1;  
    }
    return 0;
}

