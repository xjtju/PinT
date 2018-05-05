#ifndef _PinT_PFMSOLVER_Euler_H_
#define _PinT_PFMSOLVER_Euler_H_ 1

#include "pfm.h"
#include "Solver.h"
#include "PFMParams.h"

/*
 * For simple test only, directly inherited from PFMSolver, not well tuned.
 * The forward Euler is unstable for Allen-Cahn Equation
 */
class EulerSolver : public Solver {

    double beta_;       // 0.5-beta
    double dtk;         // dt*k

    void setup();
protected:
    PFMParams param;

    double *soln;   // the current solution, pointer to the grid->u_f/u_c

    size_t dsize; 
    double *soln_1;
    double *soln_2;
    double *soln_3;
    double *slns;
    
    int rfc_;
    
    void euler(long step);

public:
    EulerSolver(PinT *c, Grid *g); 
    EulerSolver(PinT *c, Grid *g, bool isFS); 

    virtual ~EulerSolver() {
        if(conf->num_std == 4) {
            free_mem(soln_1);
            free_mem(soln_2);
            free_mem(soln_3);
            free_mem(slns);
        }
    }
    virtual double* curr_solns();

    unsigned long evolve(); // Forward method iteration 

    inline void rhs() {
        if(ndim==3)      euler_rhs_ac_3d_(grid->nxyz, param.lamdaxyz, &nguard, b, soln, &dtk, &beta_);
        else if(ndim==2) euler_rhs_ac_2d_(grid->nxyz, param.lamdaxyz, &nguard, b, soln, &dtk, &beta_);
        else if(ndim==1) euler_rhs_ac_1d_(grid->nxyz, param.lamdaxyz, &nguard, b, soln, &dtk, &beta_);
    }

    inline void update() {
        if(ndim==3)      update_soln_3d_(grid->nxyz, &nguard, soln, b);
        else if(ndim==2) update_soln_2d_(grid->nxyz, &nguard, soln, b);
        else if(ndim==1) update_soln_1d_(grid->nxyz, &nguard, soln, b);
    }
};

#endif
