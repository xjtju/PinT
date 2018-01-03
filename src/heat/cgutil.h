#ifndef _CG_HEAT_UTIL_H_
#define _CG_HEAT_UTIL_H_ 1

extern "C" {

    void cg_rk2d_(int* nxyz, double* lamdaxyz, int *ng, double* r, double *x, double *b);

    void cg_xv2d_(int* nxyz, double *lamdaxyz, int *ng, double* v, double *y); 

    void cg_b2d_(int* nxyz, double* lamdaxyz, int *ng, double *x, double *b);

}
#endif
