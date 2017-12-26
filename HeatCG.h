#ifndef _PinT_HEATCG_H_
#define _PinT_HEATCG_H_ 1

#include "PBiCGStab.h"

/**
 * 1D heat equation with Crank-Nicolson
 */

class HeatCG : public PBiCGStab {

public:
    double k = 0.061644; // diffuse coefficient

    double lamda; 

    HeatCG(Grid *g, const double eps); 
    
    void cg_rk(double *r, double *x, double *b);
    void cg_Xv(double *v, double *y); 	
    void cg_b(double *x);
};

#endif
