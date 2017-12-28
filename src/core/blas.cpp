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

// s = a*v + y, v is vector, a is scalar
void blas_avpy(double* s, double a, double* v, double* y, int size){
    for(int i=0; i<size; i++){
        s[i] = a*v[i] + y[i];
    }
}

void blas_pint_sum(double *u, double *f, double *g, double *g_, double *res, int size) {
    double tmp, tmp2;
    for(int i=0; i<size; i++) {
        tmp = g[i] + f[i] - g_[i];
        tmp2 = u[i] - tmp;
        *res = *res + tmp2*tmp2;
        u[i] = tmp;
    }
}
