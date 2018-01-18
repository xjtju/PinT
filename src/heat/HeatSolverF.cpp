#include "HeatSolverF.h"
//
// the example for Fine solver and Coarse solver for heat equation.
//
// lamda = k*dt/(2*dx^2)
//  
// NOTE:
// The two classes are used only for demonstrating the usage of the framework.
// In practice, if using the same linear solver, it is not necessary to create new classes seperately for coarse and fine solver.
// 
// You can create the two instance of the base HeatSolver, 
// and setup both of them properly according to the 'isFS' parameter passed to the constructor.
//

HeatSolverF::HeatSolverF(PinT *conf, Grid *g):HeatSolver(conf,g, true) {

    lamda_x = k * conf->f_dt / (g->dx*g->dx);
    if(ndim>=2) lamda_y = k * conf->f_dt / (g->dy*g->dy);
    if(ndim>=3) lamda_z = k * conf->f_dt / (g->dz*g->dz);

    lamdaxyz[0] = lamda_x;
    lamdaxyz[1] = lamda_y;
    lamdaxyz[2] = lamda_z;

    //this->itmax = 10;  // For DEBUG
    //this->steps = 1;  // For DEBUG
    if(0 == grid->myid)
    printf("Fine : lamdax=%f, lamday=%f, lamdaz=%f \n", lamda_x, lamda_y, lamda_z);
}

HeatSolverC::HeatSolverC(PinT *conf, Grid *g):HeatSolver(conf,g, false){
    lamda_x = k * conf->c_dt / (g->dx*g->dx);
    if(ndim>=2) lamda_y = k * conf->c_dt / (g->dy*g->dy);
    if(ndim>=3) lamda_z = k * conf->c_dt / (g->dz*g->dz);

    lamdaxyz[0] = lamda_x;
    lamdaxyz[1] = lamda_y;
    lamdaxyz[2] = lamda_z;

    //this->itmax = 10;  // For DEBUG
    //this->steps = 1;  // For DEBUG
    if(0 == grid->myid)
    printf("Coarse : lamdax=%f, lamday=%f, lamdaz=%f \n", lamda_x, lamda_y, lamda_z);
}
