
#ifndef _PinT_PFMSOLVER_BD4_H_
#define _PinT_PFMSOLVER_BD4_H_ 1

#include "pfm.h"
#include "PFMSolver.h"

/*
 * The 4th order backward Euler method for Allen-Cahn Equation
 *
 * Newton-Raphson method is used to linearize the numerical formula of AC.
 */
class BD4Solver : public PFMSolver {
protected:
    void init_holder(); 

    double *soln_2;  // the holder of -F^{n-2} 
    double *soln_3;  // the holder of -F^{n-3}  
    double *soln_4;  // the holder of -F^{n-4}  

public:
    BD4Solver(PinT *c, Grid *g) : PFMSolver(c, g) {} 
    BD4Solver(PinT *c, Grid *g, bool isFS) : PFMSolver(c,g,isFS) {} 

    virtual ~BD4Solver() {
        free_mem(soln_2);
        free_mem(soln_3);
        free_mem(soln_4);
    }

    virtual void stencil() {
        if(ndim==3)      stencil_ac_3d_(grid->nxyz, lamdaxyz, &nguard, bcp, soln, &theta, &dtk, &beta_);
        else if(ndim==2) stencil_ac_2d_(grid->nxyz, lamdaxyz, &nguard, bcp, soln, &theta, &dtk, &beta_);
        else if(ndim==1) stencil_ac_1d_(grid->nxyz, lamdaxyz, &nguard, bcp, soln, &theta, &dtk, &beta_);
    }

    virtual void rhs() {
        if(ndim==3)      rhs_ac_3d_(grid->nxyz, lamdaxyz, &nguard, b, soln, soln_1, G1, &theta, &dtk, &beta_);
        else if(ndim==2) rhs_ac_2d_(grid->nxyz, lamdaxyz, &nguard, b, soln, soln_1, G1, &theta, &dtk, &beta_);
        else if(ndim==1) rhs_ac_1d_(grid->nxyz, lamdaxyz, &nguard, b, soln, soln_1, G1, &theta, &dtk, &beta_);
    }
    virtual void rhs_g1() {
        if(ndim==3)      rhs_g1_ac_3d_(grid->nxyz, lamdaxyz, &nguard, soln_1,  G1, &theta, &dtk, &beta_);
        else if(ndim==2) rhs_g1_ac_2d_(grid->nxyz, lamdaxyz, &nguard, soln_1,  G1, &theta, &dtk, &beta_);
        else if(ndim==1) rhs_g1_ac_1d_(grid->nxyz, lamdaxyz, &nguard, soln_1,  G1, &theta, &dtk, &beta_);
    }

};
#endif
