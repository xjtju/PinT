#ifndef _BLAS_CG_H_
#define _BLAS_CG_H_ 1

extern "C" {

    // 1D is not necessary to employ Fortran
    
    // 2D 
    
    // t is inner_size ,s is outer_size
    void blas_vdot_2d_( int *nxyz, int *ng, double *t, double *s, double *val); 
    //s = -av + r  or s = r - av
    void blas_avpy_2d_( int *nxyz, int *ng, double *alpha, double *s, double *r, double *v);
    // r: reverse
    void blas_avpy_2dr_(int *nxyz, int *ng, double *alpha, double *r, double *s, double *t);

    // 
    // 3D
    //
}
#endif

