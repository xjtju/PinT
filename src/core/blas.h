#ifndef PinT_BLAS_H_
#define PinT_BLAS_H_ 1

#include <math.h>
// copy s to d
void blas_cp(double* d, double* s, int size);

void blas_clear(double* d, int size);

// vector dot (v1, v2)
double blas_vdot(double* v1, double* v2, int size); 

// vector dot (v1, v2)
double blas_vdot(double* v1, double* v2, int nx, int nguard); 

// distance between two vectors 
double blas_vdist(double* v1, double* v2, int size);

// r = a*v + y, v is vector, a is scalar
void blas_avpy(double* r, double a, double* v, double* y, int size);

// r = a*v + y, v is vector, a is scalar
void blas_avpy(double* r, double a, double* v, double* y, int nx, int nguard);

void blas_pint_sum(double *u, double *f, double *g, double *g_, double *res, int size, bool fix);

void blas_pint_sum(double *u, double *f, double *g, double *g_, double *res, int nx, int nguard, bool fix);

#endif
