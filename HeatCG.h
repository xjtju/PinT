#include "PBiCGStab.h"

/**
 * 1D heat equation with Crank-Nicolson
 */

class HeatCG : public PBiCGStab {

public:
    
    double lamda = 0.061644*DT/(2*DX*DX);

    HeatCG(const int nx, const int nguard, const double eps):PBiCGStab(nx,nguard,eps){ }

    void cg_rk(double *r, double *x, double *b);
    void cg_Xv(double *v, double *y); 	
    void cg_b(double *x);
};
