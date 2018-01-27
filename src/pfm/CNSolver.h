#ifndef _PinT_PFMSOLVER_H_
#define _PinT_PFMSOLVER_H_ 1

#include "pfm.h"
#include "PBiCGStab.h"
#include "PFMParams.h"
#include "NewtonSolver.h"

/**
 * Phase Field Model using Newton-Raphson 
 *
 * Allen-Cahn Equation (AC)
 *
 * Unlike the classic Heat equation with constant diffuse coefficient, 
 * the discrete formula of the AC equation based on Crank-Nicolson is nonlinear, 
 * so Newton-Raphson method is used to linearize the numerical formula of AC, and the linear solver can be applied to AC.
 * Compared to the HEAT, it is introduced extra three variables for data structure preservation during the calculation 
 * due to its nonlinear feature.
 *
 *  The default finite difference method used in the class is Crank-Nicolson (when theta=0.5), 
 *
 * If using the default configuration, k=16000, d=1, beta=-0.128, 
 * only after 0.1 second, the system has already reached its steady state.
 * 
 */

class CNSolver : public NewtonSolver {

protected:
    PFMParams param;
    // for convenience, redefine some parameters
    double theta = 0.5;  //Crank-Nicolson; 0: Ex.Euler; 1: Im.Euler
    double beta_;       // 0.5-beta
    double dtk;         // dt*k

    void setup();

    virtual void stencil() {
        if(ndim==3)      stencil_ac_3d_(grid->nxyz, param.lamdaxyz, &nguard, A, soln, &theta, &dtk, &beta_);
        else if(ndim==2) stencil_ac_2d_(grid->nxyz, param.lamdaxyz, &nguard, A, soln, &theta, &dtk, &beta_);
        else if(ndim==1) stencil_ac_1d_(grid->nxyz, param.lamdaxyz, &nguard, A, soln, &theta, &dtk, &beta_);
    }

    virtual void rhs() {
        if(ndim==3)      rhs_ac_3d_(grid->nxyz, param.lamdaxyz, &nguard, b, soln, soln_1, G1, &theta, &dtk, &beta_);
        else if(ndim==2) rhs_ac_2d_(grid->nxyz, param.lamdaxyz, &nguard, b, soln, soln_1, G1, &theta, &dtk, &beta_);
        else if(ndim==1) rhs_ac_1d_(grid->nxyz, param.lamdaxyz, &nguard, b, soln, soln_1, G1, &theta, &dtk, &beta_);
    }
    virtual void rhs_g1() {
        if(ndim==3)      rhs_g1_ac_3d_(grid->nxyz, param.lamdaxyz, &nguard, soln_1,  G1, &theta, &dtk, &beta_);
        else if(ndim==2) rhs_g1_ac_2d_(grid->nxyz, param.lamdaxyz, &nguard, soln_1,  G1, &theta, &dtk, &beta_);
        else if(ndim==1) rhs_g1_ac_1d_(grid->nxyz, param.lamdaxyz, &nguard, soln_1,  G1, &theta, &dtk, &beta_);
    }

public:

    CNSolver(PinT *c, Grid *g); 
    CNSolver(PinT *c, Grid *g, bool isFS); 

    virtual ~CNSolver() {
    }
};

#endif
