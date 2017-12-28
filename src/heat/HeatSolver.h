#ifndef _PinT_HEATSOLVER_H_
#define _PinT_HEATSOLVER_H_ 1

#include "PBiCGStab.h"

/**
 * 1D heat equation with Crank-Nicolson
 */

class HeatSolver : public PBiCGStab {

public:
    double k = 0.061644; // diffuse coefficient

    double lamda; 

    HeatSolver(PinT *c, Grid *g); 
    
    void cg_rk(double *r, double *x, double *b);
    void cg_Xv(double *v, double *y); 	
    void cg_b(double *x);
};

#endif
