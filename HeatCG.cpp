#include "HeatCG.h"

HeatCG::HeatCG(Grid *g, const double eps) : PBiCGStab(g,eps){
    lamda = k * g->dt / (2*g->dx*g->dx); 
}

//calcaluate the residual r = b - Ax
void HeatCG::cg_rk(double *r, double *x, double *b){
    for(int i=nguard; i<nx+nguard; i++){
        double ax = -lamda*x[i-1] + (1+2*lamda)*x[i] - lamda*x[i+1];  
        r[i] = b[i] - ax; 
    }
}

// matrix * vector,  the stencil matrix 
// v = Ay
void HeatCG::cg_Xv(double* v, double *y) {
    for(int i=nguard; i<nx+nguard; i++){
        v[i] = -lamda*y[i-1] + (1+2*lamda)*y[i] - lamda*y[i+1];
    }
}

void HeatCG::cg_b(double *x){
    for(int i=nguard; i<nx+nguard; i++){
        b[i] = lamda*x[i-1] + (1-2*lamda)*x[i] + lamda*x[i+1]; 
    }
}
