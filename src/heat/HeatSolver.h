#ifndef _PinT_HEATSOLVER_H_
#define _PinT_HEATSOLVER_H_ 1

#include "PBiCGStab.h"
#include "cgutil.h"
/**
 * 1D/2D heat equation with Crank-Nicolson
 */

class HeatSolver : public PBiCGStab {

protected:

    void setup();
public:
    double k = 0.061644; // diffuse coefficient

    // lamda = k*dt/dx^2   
    // when dx, dy, dz are not equal with each other, the lamdas are unequall neither.
    double lamda_x;
    double lamda_y;
    double lamda_z;
    double &lamda = lamda_x;

    double lamdaxyz[3];

    HeatSolver(PinT *c, Grid *g); 
    HeatSolver(PinT *c, Grid *g, bool isFS); 
    
    void cg_rk1d(double *r, double *x, double *b);
    void cg_rk2d(double *r, double *x, double *b);

    void cg_Xv1d(double *v, double *y); 	
    void cg_Xv2d(double *v, double *y); 	

    void cg_b1d(double *x);
    void cg_b2d(double *x);

    inline void sor2_core_1d(double *p_, double *p, int *color){
         sor2_core_1d_(grid->nxyz, lamdaxyz, &nguard, p_, p, color, &sor_omg) ;
    }

    // SOR for  Ap_=p
    inline void sor2_core_2d(double *p_, double *p, int *color) {
         sor2_core_2d_(grid->nxyz, lamdaxyz, &nguard, p_, p, color, &sor_omg) ;
    }
};

#endif
