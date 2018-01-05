#ifndef _PinT_BLAS_CG_H_
#define _PinT_BLAS_CG_H_ 1

extern "C" {

    // 1D,  
    
    void blas_avpy_1d_( int *nxyz, int *ng, double *alpha, double *s, double *r, double *v);

    void cg_xi_1d_( int *nxyz, int *ng, double *x, double *y, double *z, double *alpha, double *omega);
    
    void cg_direct_1d_( int *nxyz, int *ng, double* p, double* r, double* v, double* beta, double* omega);

    // 2D 
    // 

    //s = -av + r  or s = r - av
    void blas_avpy_2d_( int *nxyz, int *ng, double *alpha, double *s, double *r, double *v);

    // x = x + ay + bz  
    void cg_xi_2d_( int *nxyz, int *ng, double *x, double *y, double *z, double *alpha, double *omega);
    
    //p = r + beta * ( p - omg * v )
    void cg_direct_2d_( int *nxyz, int *ng, double* p, double* r, double* v, double* beta, double* omega);

    // 
    // 3D
    //
}
#endif

