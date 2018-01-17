#ifndef _PinT_PFMSOLVER_H_
#define _PinT_PFMSOLVER_H_ 1

#include "PBiCGStab.h"
#include "pfm.h"

/**
 * Phase Field Model using Crank-Nicolson
 */

class PFMSolver : public Solver {

protected:
    void setup(PinT *c, Grid *g);

public:

    double k;     // interfacial width related
    double d;     // diffuse coefficient 
    double beta;

    double theta = 0.5;  //Crank-Nicolson
    int newton_itmax = 5;

    double lamda_x;

    double *F;   // the current solution, pointer to the grid->u_f/u_c
    double *F_;  // the holder of -F^{k-1} in Newton's method when applying to nonlinear systems of equations 
    double *G1;  // the partial RHS of Newton's method for PFM  
    double *unk; // the unknown 'x' of Ax=b, that is (Xn+1 - Xn) for Newton's method 
    double *b;   //RHS
    double *bcp; //stencil matrix

    double A0, A1_, A1; // Jacobi matrix
    double dtk;         // dt*k
    double beta_;       // 0.5-beta

    double ls_eps;
    double ls_itmax;

    PBiCGStab *hypre;  // linear solver 

    PFMSolver(PinT *c, Grid *g); 
    PFMSolver(PinT *c, Grid *g, bool isFS); 

    ~PFMSolver() {
        free_mem(F_);
        free_mem(G1);
        free_mem(unk);
        free_mem(b);
        free_mem(bcp);
    }

    void evolve();   // evolve for New-Raphson method
    void newton_raphson(); 
   
    inline void stencil() {
        if(ndim==1) stencil_ac_1d_(grid->nxyz, &lamda_x, &nguard, bcp, F, &theta, &dtk, &beta_);
    }
    void rhs();
};

int pfm_inih(void* obj, const char* section, const char* name, const char* value);

#endif
