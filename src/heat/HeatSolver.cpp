#include "HeatSolver.h"

HeatSolver::HeatSolver(PinT *c, Grid *g) : PBiCGStab(c,g){
    lamda_x = k * conf->f_dt / (2*g->dx*g->dx);
    lamda_y = k * conf->f_dt / (2*g->dy*g->dy);
    lamda_z = k * conf->f_dt / (2*g->dz*g->dz);

    lamdaxyz[0] = lamda_x;
    lamdaxyz[1] = lamda_y;
    lamdaxyz[2] = lamda_z;

    this->eps = 1.0e-6;
    this->itmax = 10;
}

//calcaluate the residual r = b - Ax
void HeatSolver::cg_rk1d(double *r, double *x, double *b){
    for(int i=nguard; i<nx+nguard; i++){
        double ax = -lamda*x[i-1] + (1+2*lamda)*x[i] - lamda*x[i+1];  
        r[i] = b[i] - ax;
    }
}

// matrix * vector,  the stencil matrix 
// v = Ay
void HeatSolver::cg_Xv1d(double* v, double *y) {
    for(int i=nguard; i<nx+nguard; i++){
        v[i] = -lamda*y[i-1] + (1+2*lamda)*y[i] - lamda*y[i+1];
    }
}

void HeatSolver::cg_b1d(double *x){
    for(int i=nguard; i<nx+nguard; i++){
        b[i] = lamda*x[i-1] + (1-2*lamda)*x[i] + lamda*x[i+1]; 
    }
}

//
// 2D
//

void HeatSolver::cg_rk2d(double *r, double *x, double *b){

    cg_rk2d_(grid->nxyz, lamdaxyz, &nguard, r, x, b);
}

// matrix * vector,  the stencil matrix 
// v = Ay
void HeatSolver::cg_Xv2d(double* v, double *y) {

    cg_xv2d_(grid->nxyz, lamdaxyz, &nguard, v, y);
}

void HeatSolver::cg_b2d(double *x){
    cg_b2d_(grid->nxyz, lamdaxyz, &nguard, x, b);
}
