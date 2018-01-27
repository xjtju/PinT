#ifndef _PinT_PFMPARAMS_H_
#define _PinT_PFMPARAMS_H_ 1

#include "PinT.h"
#include "Grid.h"

//problem specific parameter parsing
int pfm_inih(void* obj, const char* section, const char* name, const char* value);

struct PFMParams {

    double k;     // interfacial width related
    double d;     // diffuse coefficient 
    double beta;  // potential engrgy parameter
    double beta_;       // 0.5-beta

    double theta = 0.5;  //Crank-Nicolson; 0: Ex.Euler; 1: Im.Euler

    double dtk;         // dt*k
    double lamda_x;     // d*dt/(dx**2)
    double lamda_y = 0.0;     
    double lamda_z = 0.0;    
    double lamdaxyz[3];

    int newton_itmax = 10;
    
    int ndim; 

    int fsolver = 0;
    int csolver = 0;
    
    void init(PinT* conf, Grid *g){
        conf->init_module(this, pfm_inih);
        beta_ = 0.5 - beta;
        ndim = conf->ndim;
    }

    void init(PinT* conf, Grid *g, bool isFS){
        init(conf, g);

        if(isFS) {
            lamda_x = d*conf->f_dt/(g->dx*g->dx);   
            if(ndim>=2) lamda_y = d * conf->f_dt / (g->dy*g->dy);
            if(ndim>=3) lamda_z = d * conf->f_dt / (g->dz*g->dz);
            dtk = conf->f_dt*k;
        }else {
            lamda_x = d*conf->c_dt/(g->dx*g->dx);   
            if(ndim>=2) lamda_y = d * conf->c_dt / (g->dy*g->dy);
            if(ndim>=3) lamda_z = d * conf->c_dt / (g->dz*g->dz);
            dtk = conf->c_dt*k; 
        }

        lamdaxyz[0] = lamda_x;
        lamdaxyz[1] = lamda_y;
        lamdaxyz[2] = lamda_z;
    }

    void print() {
        printf("PFM init parameter : \n");  
        printf("        D    = %f \n", d);
        printf("        k    = %f \n", k);
        printf("        beta = %f \n", beta);
        printf("  newton_iter= %d \n", newton_itmax);
        printf("        theta= %f , (0.0: Ex.Euler; 0.5: Crank-Nicolson; 1.0: Im.Euler.)\n\n", theta);
        printf("  Fine Solver= %s \n", getSolvername(fsolver));
        printf("  Coar Solver= %s \n", getSolvername(csolver));
    
    }
    void printLamda(const char* flag) {
        printf("%s : dtk=%f, lamda_x=%f, lamday=%f, lamdaz=%f \n", flag, dtk, lamda_x, lamda_y, lamda_z );
    }

    const char* CN  = "CN";
    const char* BD  = "BD4";
    const char* EU  = "Euler";
    const int ID_CN = 0;
    const int ID_BD = 1;
    const int ID_EU = 2;

    const char* getSolvername(int s) {
        if(s==ID_CN) return CN;
        if(s==ID_BD) return BD;
        if(s==ID_EU) return EU;
    }
};

#endif

