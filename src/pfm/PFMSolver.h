#ifndef _PinT_PFMSOLVER_H_
#define _PinT_PFMSOLVER_H_ 1

#include "PBiCGStab.h"

/**
 * Phase Field Model using Crank-Nicolson
 */

class PFMSolver : public PBiCGStab {

protected:
    void setup();

public:

    double k;     // interfacial width related
    double d;     // diffuse coefficient 
    double beta;

    double theta = 0.5;  //Crank-Nicolson
    int newton_itmax = 5;

    double lamda_x;

    double *F_;  // the holder of -F^{k-1} in Newton's method when applying to nonlinear systems of equations 
    double *F;
    double *G1;  // the pointer to the starting value 
    double *unk; //Xn+1 - Xn

    double A0, A1_, A1; // Jacobi matrix
    double dtk;         // dt*k
    double beta_;       // 0.5-beta

    PFMSolver(PinT *c, Grid *g); 
    PFMSolver(PinT *c, Grid *g, bool isFS); 
     
    void newton_raphson(); 
    void update(); 
    double* fetch(); 
    void prepare();
    void evolve();   // overwrite the default evolve for New-Raphson method

    void cg_rk1d(double *r, double *x, double *b);
    void cg_rk2d(double *r, double *x, double *b);
    void cg_rk3d(double *r, double *x, double *b);

    void cg_Xv1d(double *v, double *y); 	
    void cg_Xv2d(double *v, double *y); 	
    void cg_Xv3d(double *v, double *y); 	

    void cg_b1d(double *x);
    void cg_b2d(double *x);
    void cg_b3d(double *x);
    
   
    // the preconditioner is not used at current
    // SOR for  Ap_=p
    inline void sor2_core_1d(double *p_, double *p, int *color){ }
    inline void sor2_core_2d(double *p_, double *p, int *color) { }
    
};

int pfm_inih(void* obj, const char* section, const char* name, const char* value);

#endif
