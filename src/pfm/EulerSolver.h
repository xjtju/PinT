#ifndef _PinT_PFMSOLVER_Euler_H_
#define _PinT_PFMSOLVER_Euler_H_ 1

#include "pfm.h"
#include "PFMSolver.h"

/*
 * For simple test only, directly inherited from PFMSolver, not well tuned.
 * The forward Euler is unstable for Allen-Cahn Equation
 */
class EulerSolver : public PFMSolver {
protected:
    void euler();

public:
    EulerSolver(PinT *c, Grid *g); 
    EulerSolver(PinT *c, Grid *g, bool isFS); 

    virtual ~EulerSolver() {}

    void evolve(); // Forward method iteration 

    inline void rhs() {
        if(ndim==3)      euler_rhs_ac_3d_(grid->nxyz, lamdaxyz, &nguard, b, soln, &theta, &dtk, &beta_);
        else if(ndim==2) euler_rhs_ac_2d_(grid->nxyz, lamdaxyz, &nguard, b, soln, &theta, &dtk, &beta_);
        else if(ndim==1) euler_rhs_ac_1d_(grid->nxyz, lamdaxyz, &nguard, b, soln, &theta, &dtk, &beta_);
    }

    inline void update() {
        if(ndim==3)      update_ac_3d_(grid->nxyz, &nguard, soln, b);
        else if(ndim==2) update_ac_2d_(grid->nxyz, &nguard, soln, b);
        else if(ndim==1) update_ac_1d_(grid->nxyz, &nguard, soln, b);
    }
};

#endif
