#include "HeatSolverF.h"
//
// the example for Fine solver and Coarse solver for heat equation.
//
// lamda = k*dt/(2*dx^2)
//

HeatSolverF::HeatSolverF(PinT *conf, Grid *g):HeatSolver(conf,g, true) {

    lamda_x = k * conf->f_dt / (2*g->dx*g->dx);
    if(ndim>=2) lamda_y = k * conf->f_dt / (2*g->dy*g->dy);
    if(ndim>=3) lamda_z = k * conf->f_dt / (2*g->dz*g->dz);

    lamdaxyz[0] = lamda_x;
    lamdaxyz[1] = lamda_y;
    lamdaxyz[2] = lamda_z;

    //this->itmax = 10;  // For DEBUG
    //this->steps = 1;  // For DEBUG
    if(0 == grid->myid)
    printf("Fine : lamdax=%f, lamday=%f, lamdaz=%f \n", lamda_x, lamda_y, lamda_z);
}

HeatSolverC::HeatSolverC(PinT *conf, Grid *g):HeatSolver(conf,g, false){
    lamda_x = k * conf->c_dt / (2*g->dx*g->dx);
    if(ndim>=2) lamda_y = k * conf->c_dt / (2*g->dy*g->dy);
    if(ndim>=3) lamda_z = k * conf->c_dt / (2*g->dz*g->dz);

    lamdaxyz[0] = lamda_x;
    lamdaxyz[1] = lamda_y;
    lamdaxyz[2] = lamda_z;

    //this->itmax = 10;  // For DEBUG
    //this->steps = 1;  // For DEBUG
    if(0 == grid->myid)
    printf("Coarse : lamdax=%f, lamday=%f, lamdaz=%f \n", lamda_x, lamda_y, lamda_z);
}
