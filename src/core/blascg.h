#ifndef _BLAS_CG_H_
#define _BLAS_CG_H_ 1

extern "C" {

    // 1D is not necessary to employ Fortran
    
    // 2D 
    
    // t is inner_size ,s is outer_size
    void blas_vdot_2d_( int *nxyz, int *ng, double *t, double *s, double *val); 

    //s = -av + r  or s = r - av
    void blas_avpy_2d_( int *nxyz, int *ng, double *alpha, double *s, double *r, double *v);

    // x = x + ay + bz  
    void cg_xi2d_( int *nxyz, int *ng, double *x, double *y, double *z, double *alpha, double *omega);
    
    //p = r + beta * ( p - omg * v )
    void cg_direct2d_( int *nxyz, int *ng, double* p, double* r, double* v, double* beta, double* omega);

    // 
    // 3D
    //
}
#endif

