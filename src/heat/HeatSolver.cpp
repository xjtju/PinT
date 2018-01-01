#include "HeatSolver.h"

HeatSolver::HeatSolver(PinT *c, Grid *g) : PBiCGStab(c,g){ }

//calcaluate the residual r = b - Ax
void HeatSolver::cg_rk(double *r, double *x, double *b){
    for(int i=nguard; i<nx+nguard; i++){
        double ax = -lamda*x[i-1] + (1+2*lamda)*x[i] - lamda*x[i+1];  
        r[i] = b[i] - ax;
        //printf("%f ",r[i]);
    }
}

// matrix * vector,  the stencil matrix 
// v = Ay
void HeatSolver::cg_Xv(double* v, double *y) {
    for(int i=nguard; i<nx+nguard; i++){
        v[i] = -lamda*y[i-1] + (1+2*lamda)*y[i] - lamda*y[i+1];
    }
}

void HeatSolver::cg_b(double *x){
    for(int i=nguard; i<nx+nguard; i++){
        b[i] = lamda*x[i-1] + (1-2*lamda)*x[i] + lamda*x[i+1]; 
    }
}
