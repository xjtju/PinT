#include "HeatSolverF.h"


HeatSolverF::HeatSolverF(PinT *conf, Grid *g):HeatSolver(conf,g) {

    lamda_x = k * conf->f_dt / (2*g->dx*g->dx);
    lamda_y = k * conf->f_dt / (2*g->dy*g->dy);
    lamda_z = k * conf->f_dt / (2*g->dz*g->dz);

    lamdaxyz[0] = lamda_x;
    lamdaxyz[1] = lamda_y;
    lamdaxyz[2] = lamda_z;

    this->steps = conf->f_steps;  

    printf("fine   : lamdax=%f, lamday=%f \n", lamda_x, lamda_y);
}

double* HeatSolverF::fetch(){
     return grid->u_f;
}

HeatSolverC::HeatSolverC(PinT *conf, Grid *g):HeatSolver(conf,g){
    lamda_x = k * conf->c_dt / (2*g->dx*g->dx);
    lamda_y = k * conf->c_dt / (2*g->dy*g->dy);
    lamda_z = k * conf->c_dt / (2*g->dz*g->dz);

    lamdaxyz[0] = lamda_x;
    lamdaxyz[1] = lamda_y;
    lamdaxyz[2] = lamda_z;

    this->steps = conf->c_steps; 

    printf("Coarse : lamdax=%f, lamday=%f \n", lamda_x, lamda_y);
}
double* HeatSolverC::fetch(){
     return grid->u_c;
}
