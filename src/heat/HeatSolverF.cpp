#include "HeatSolverF.h"
//
// the example for Fine solver and Coarse solver for heat equation.
//

HeatSolverF::HeatSolverF(PinT *conf, Grid *g):HeatSolver(conf,g, true) {

    lamda_x = k * conf->f_dt / (2*g->dx*g->dx);
    lamda_y = k * conf->f_dt / (2*g->dy*g->dy);
    lamda_z = k * conf->f_dt / (2*g->dz*g->dz);

    lamdaxyz[0] = lamda_x;
    lamdaxyz[1] = lamda_y;
    lamdaxyz[2] = lamda_z;

    //this->itmax = 10;  // For DEBUG
    //this->steps = 10;  // For DEBUG
    if(0 == grid->myid)
    printf("fine   : lamdax=%f, lamday=%f \n", lamda_x, lamda_y);
}

HeatSolverC::HeatSolverC(PinT *conf, Grid *g):HeatSolver(conf,g, false){
    lamda_x = k * conf->c_dt / (2*g->dx*g->dx);
    lamda_y = k * conf->c_dt / (2*g->dy*g->dy);
    lamda_z = k * conf->c_dt / (2*g->dz*g->dz);

    lamdaxyz[0] = lamda_x;
    lamdaxyz[1] = lamda_y;
    lamdaxyz[2] = lamda_z;

    //this->itmax = 10;  // For DEBUG
    //this->steps = 100;  // For DEBUG
    if(0 == grid->myid)
    printf("Coarse : lamdax=%f, lamday=%f \n", lamda_x, lamda_y);
}
