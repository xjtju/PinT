#ifndef _FORTRAN_SOR_H_
#define _FORTRAN_SOR_H_ 1

extern "C" {

    void sor2_core_1d_(int* nxyz, int *ng, double *p_, double *p, double *A, int *color, double *omg, double *err) ; 
    void sor2_core_2d_(int* nxyz, int *ng, double *p_, double *p, double *A, int *color, double *omg, double *err) ;
    void sor2_core_3d_(int* nxyz, int *ng, double *p_, double *p, double *A, int *color, double *omg, double *err) ;
}
#endif
