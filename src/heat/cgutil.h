#ifndef _CG_HEAT_UTIL_H_
#define _CG_HEAT_UTIL_H_ 1

//
//  HEAT stencil by Fortran
//
extern "C" {

    void cg_rk2d_(int* nxyz, double* lamdaxyz, int *ng, double* r, double *x, double *b);
    void cg_xv2d_(int* nxyz, double *lamdaxyz, int *ng, double* v, double *y); 
    void cg_b2d_(int* nxyz, double* lamdaxyz, int *ng, double *x, double *b);

    void cg_rk3d_(int* nxyz, double* lamdaxyz, int *ng, double* r, double *x, double *b);
    void cg_xv3d_(int* nxyz, double *lamdaxyz, int *ng, double* v, double *y); 
    void cg_b3d_(int* nxyz, double* lamdaxyz, int *ng, double *x, double *b);

    void sor2_core_1d_(int* nxyz, double* lamdaxyz, int *ng, double *p_, double *p, int *color, double *omg) ;
    void sor2_core_2d_(int* nxyz, double* lamdaxyz, int *ng, double *p_, double *p, int *color, double *omg) ;
}
#endif
