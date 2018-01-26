
#ifndef _PinT_PFMSOLVER_BD4_H_
#define _PinT_PFMSOLVER_BD4_H_ 1

#include "pfm.h"
#include "PFMSolver.h"

/*
 * The 4th order backward differentiation formula for Allen-Cahn Equation
 * 
 * Because the time integration part is inherited from PFMSolver, 
 * just only BD4-specific calculations have to be implemented.     
 */
class BD4Solver : public PFMSolver {
protected:
    void create_holder(); 

    double *soln_2;  // the holder of -F^{n-2} 
    double *soln_3;  // the holder of -F^{n-3}  
    double *soln_4;  // the holder of -F^{n-4}  

public:
    BD4Solver(PinT *c, Grid *g) : PFMSolver(c, g) {
        create_holder();
    } 
    BD4Solver(PinT *c, Grid *g, bool isFS) : PFMSolver(c,g,isFS) {
        create_holder();
    } 

    virtual ~BD4Solver() {
        free_mem(soln_2);
        free_mem(soln_3);
        free_mem(soln_4);
        if(grid->myid==0 && conf->verbose)
            printf("INFO: the memory allocate by BD4Solver has been released.\n");
    }

    // It is not necessary to make the following functions 'virtual', 
    // just for clearly marking that the default implementation of the superclass has been overwritten.
    virtual void update_holder();
    virtual void init_holder();

    virtual void rhs() {
        if(ndim==3)      rhs_ac_bd4_3d_(grid->nxyz, lamdaxyz, &nguard, b, soln, soln_1, G1, &dtk, &beta_);
        else if(ndim==2) rhs_ac_bd4_2d_(grid->nxyz, lamdaxyz, &nguard, b, soln, soln_1, G1, &dtk, &beta_);
        else if(ndim==1) rhs_ac_bd4_1d_(grid->nxyz, lamdaxyz, &nguard, b, soln, soln_1, G1, &dtk, &beta_);
    }

    virtual void stencil() {
        if(ndim==3)      stencil_ac_bd4_3d_(grid->nxyz, lamdaxyz, &nguard, bcp, soln, &dtk, &beta_);
        else if(ndim==2) stencil_ac_bd4_2d_(grid->nxyz, lamdaxyz, &nguard, bcp, soln, &dtk, &beta_);
        else if(ndim==1) stencil_ac_bd4_1d_(grid->nxyz, lamdaxyz, &nguard, bcp, soln, &dtk, &beta_);
    }

    virtual void rhs_g1() {
        if(ndim==3)      rhs_g1_ac_bd4_3d_(grid->nxyz, &nguard, soln_1, soln_2, soln_3, soln_4, G1, &dtk, &beta_);
        else if(ndim==2) rhs_g1_ac_bd4_2d_(grid->nxyz, &nguard, soln_1, soln_2, soln_3, soln_4, G1, &dtk, &beta_);
        else if(ndim==1) rhs_g1_ac_bd4_1d_(grid->nxyz, &nguard, soln_1, soln_2, soln_3, soln_4, G1, &dtk, &beta_);
    }

};
#endif
