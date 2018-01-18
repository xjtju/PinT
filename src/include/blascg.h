#ifndef _PinT_BLAS_CG_H_
#define _PinT_BLAS_CG_H_ 1

extern "C" {

    // 1D,  
    
    void blas_avpy_1d_( int *nxyz, int *ng, double *alpha, double *s, double *r, double *v);

    void cg_xi_1d_( int *nxyz, int *ng, double *x, double *y, double *z, double *alpha, double *omega);
    
    void cg_direct_1d_( int *nxyz, int *ng, double* p, double* r, double* v, double* beta, double* omega);

    //calcaluate the residual r = b - Ax
    void cg_rk1d_( int *nxyz, int *ng, double *r, double *x, double *b, double *bcp);
    // v = Ax
    void cg_ax1d_( int *nxyz, int *ng, double *v, double *x, double *bcp);

    // 2D 
    // 

    //s = -av + r  or s = r - av
    void blas_avpy_2d_( int *nxyz, int *ng, double *alpha, double *s, double *r, double *v);

    // x = x + ay + bz  
    void cg_xi_2d_( int *nxyz, int *ng, double *x, double *y, double *z, double *alpha, double *omega);
    
    //p = r + beta * ( p - omg * v )
    void cg_direct_2d_( int *nxyz, int *ng, double* p, double* r, double* v, double* beta, double* omega);

    void cg_rk2d_( int *nxyz, int *ng, double *r, double *x, double *b, double *bcp);
    void cg_ax2d_( int *nxyz, int *ng, double *v, double *x, double *bcp);

    // 
    // 3D
    //
    void blas_avpy_3d_( int *nxyz, int *ng, double *alpha, double *s, double *r, double *v);
    void cg_xi_3d_( int *nxyz, int *ng, double *x, double *y, double *z, double *alpha, double *omega);
    void cg_direct_3d_( int *nxyz, int *ng, double* p, double* r, double* v, double* beta, double* omega);

    void cg_rk3d_( int *nxyz, int *ng, double *r, double *x, double *b, double *bcp);
    void cg_ax3d_( int *nxyz, int *ng, double *v, double *x, double *bcp);

}
#endif

