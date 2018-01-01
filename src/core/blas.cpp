#include "blas.h"


//copy vector s to d
void blas_cp(double* d, double* s, int size){
    for(int i=0; i<size; i++) {
        d[i] = s[i];
    }
}

//vector dot (v1, v2)
double blas_vdot(double* v1, double* v2, int size){
    double r = 0;
    for(int i=0; i<size; i++) {
        r = r + v1[i]*v2[i];
    }
    return r;
}

//vector dot (v1, v2)
double blas_vdot(double* v1, double* v2, int nx, int nguard){
    double r = 0;
    for(int i=nguard; i<nx+nguard; i++) {
        r = r + v1[i]*v2[i];
    }
    return r;
}

double blas_vdist(double* v1, double* v2, int size){
    double dist = 0.0;
    for(int i=0; i<size; i++){
        dist = dist + (v1[i]-v2[i])*(v1[i]-v2[i]); 
    }
    return sqrt(dist); 
}

// s = a*v + y, v is vector, a is scalar
void blas_avpy(double* s, double a, double* v, double* y, int size){
    for(int i=0; i<size; i++){
        s[i] = a*v[i] + y[i];
    }
}

// s = a*v + y, v is vector, a is scalar
void blas_avpy(double* s, double a, double* v, double* y, int nx, int nguard){
    for(int i=nguard; i<nx+nguard; i++){
        s[i] = a*v[i] + y[i];
    }
}

/**
 * Due to the "Machine Epsilon" or rounding error, residual calculation is very important.
 * Sometimes, in theory, the residual should be ZERO, but in practice, the calculation value is not ZERO, 
 * despite it is very small, it will has an unignorable impact on convergency due to  accumulating effect.  
 */
void blas_pint_sum(double *u, double *f, double *g, double *g_, double *res, int size, bool fix)  {
    double tmp = 0.0;
    double tmp2 = 0.0;
    for(int i=0; i<size; i++) {
        tmp  = f[i] + (g[i] - g_[i]);
        tmp2 = u[i] - tmp;
        *res = *res + tmp2*tmp2;
        u[i] = tmp;
    }
    if(fix) {
        if (*res < 1.0e-12) *res = 0.0;
    }
}

void blas_pint_sum(double *u, double *f, double *g, double *g_, double *res, int nx, int nguard, bool fix)  {
    double tmp = 0.0;
    double tmp2 = 0.0;
    for(int i=nguard; i<nx+nguard; i++) {
        tmp  = f[i] + (g[i] - g_[i]);
        tmp2 = u[i] - tmp;
        *res = *res + tmp2*tmp2;
        u[i] = tmp;
    }
    if(fix) {
        if (*res < 1.0e-12) *res = 0.0;
    }
}
