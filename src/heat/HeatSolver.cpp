#include "HeatSolver.h"

HeatSolver::HeatSolver(PinT *c, Grid *g) : PBiCGStab(c,g){
    setup();
}

HeatSolver::HeatSolver(PinT *c, Grid *g, bool isFS) : PBiCGStab(c,g,isFS){
    setup();
}

// set diffuse coefficient and tune the default parameter, problem specific
void HeatSolver::setup(){
    this->eps = 1.0e-6;

    if(ndim==1) k = 0.061644; 
    if(ndim==2) k = 0.061644;
}

//
// 1D, not used Fortran
//

//calcaluate the residual r = b - Ax
void HeatSolver::cg_rk1d(double *r, double *x, double *b){
    int ind;
    for(int i=nguard; i<nx+nguard; i++){
        ind = i; //grid->getInnerIdx(i);
        double ax = -lamda*x[i-1] + (1+2*lamda)*x[i] - lamda*x[i+1];  
        r[ind] = b[ind] - ax;
    }
}

// matrix * vector,  the stencil matrix 
// v = Ay
void HeatSolver::cg_Xv1d(double* v, double *y) {
    int ind;
    for(int i=nguard; i<nx+nguard; i++){
        ind = i; //grid->getInnerIdx(i);
        v[ind] = -lamda*y[i-1] + (1+2*lamda)*y[i] - lamda*y[i+1];
    }
}

void HeatSolver::cg_b1d(double *x){
    int ind;
    for(int i=nguard; i<nx+nguard; i++){
        ind = i; //grid->getInnerIdx(i);
        b[ind] = lamda*x[i-1] + (1-2*lamda)*x[i] + lamda*x[i+1]; 
    }
}

//
// 2D, used Fortran
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


