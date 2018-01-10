#ifndef PinT_BLAS_H_
#define PinT_BLAS_H_ 1

extern "C" {

    void blas_clear_(double* d, size_t *size);
    void blas_cp_(double *d, double* s, size_t *size);

    // vector dot (v1, v2)
    void blas_dot_1d_( int *nxyz, int *ng, double *t, double *s, double *val); 
    void blas_dot_2d_( int *nxyz, int *ng, double *t, double *s, double *val); 
    void blas_dot_3d_( int *nxyz, int *ng, double *t, double *s, double *val); 

    void blas_vdist_1d_( int *nxyz, int *ng, double* v1, double* v2, double *val);
    void blas_vdist_2d_( int *nxyz, int *ng, double* v1, double* v2, double *val);
    void blas_vdist_3d_( int *nxyz, int *ng, double* v1, double* v2, double *val);

    /**
     * Due to the "Machine Epsilon" or rounding error, residual calculation is very important.
     * Sometimes, in theory, the residual should be ZERO, but in practice, the calculation value is not ZERO, 
     * despite it is very small, it will has an unignorable impact on convergency due to  accumulating effect.  
     */
     void blas_pint_sum_1d_(int *nxyz, int *ng, double *u, double *f, double *g, double *g_, double *res, double *sml);  
     void blas_pint_sum_2d_(int *nxyz, int *ng, double *u, double *f, double *g, double *g_, double *res, double *sml);  
     void blas_pint_sum_3d_(int *nxyz, int *ng, double *u, double *f, double *g, double *g_, double *res, double *sml);  
}
#endif
