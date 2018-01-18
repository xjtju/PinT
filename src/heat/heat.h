#ifndef _CG_HEAT_UTIL_H_
#define _CG_HEAT_UTIL_H_ 1

//
//  HEAT stencil by Fortran
//
extern "C" {

    void sor2_core_1d_(int* nxyz, double* lamdaxyz, int *ng, double *p_, double *p, int *color, double *omg) ; // not used
    void sor2_core_2d_(int* nxyz, double* lamdaxyz, int *ng, double *p_, double *p, int *color, double *omg) ; // not used

    void rhs_heat_1d_(    int *nxyz, double *lamda, int *ng, double *soln, double *b);
    void stencil_heat_1d_(int *nxyz, double *lamda, int *ng, double *soln, double *bcp);

    void rhs_heat_2d_(    int *nxyz, double *lamda, int *ng, double *soln, double *b);
    void stencil_heat_2d_(int *nxyz, double *lamda, int *ng, double *soln, double *bcp);

    void rhs_heat_3d_(    int *nxyz, double *lamda, int *ng, double *soln, double *b);
    void stencil_heat_3d_(int *nxyz, double *lamda, int *ng, double *soln, double *bcp);
}
#endif
