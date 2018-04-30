
#ifndef _PinT_PFMSOLVER_BD4_H_
#define _PinT_PFMSOLVER_BD4_H_ 1

#include "pfm.h"
#include "PFMParams.h"
#include "NewtonSolver.h"
#include "Driver.h"
/*
 * The 4th order backward differentiation formula for Allen-Cahn Equation
 * 
 * Because the time integration part is inherited from PFMSolver, 
 * just only BD4-specific calculations have to be implemented.     
 */
class BD4Solver : public NewtonSolver {
private:
    void setup();
    void create_holder(); 

    double *soln_2;  // the holder of U^{l-2} 
    double *soln_3;  // the holder of U^{l-3}  
    double *soln_4;  // the holder of U^{l-4}, where 'l' is the index of timestep  

    double *prevslns;
    double *slns;
    double *sbuf;
    double *rbuf;
    
    PFMParams param;
    double dtk;
    double beta_;
    size_t dsize;

public:
    BD4Solver(PinT *c, Grid *g) ; 
    BD4Solver(PinT *c, Grid *g, bool isFS) ; 

    virtual ~BD4Solver() {
        free_mem(soln_2);
        free_mem(soln_3);
        free_mem(soln_4);
        free_mem(sbuf);
        free_mem(rbuf);
        free_mem(prevslns);
        free_mem(slns);
        if(grid->myid==0 && conf->verbose)
            printf("INFO: the memory allocate by BD4Solver has been released.\n");
    }

protected:
    // It is not necessary to make the following functions 'virtual', 
    // just for clearly marking that the default implementation of the superclass has been overwritten.
    virtual void update_holder();
    virtual void init_holder();


     virtual double* prev_solns() {
         return prevslns;
     }

     virtual double* curr_solns();

     // preserve the solution of coarse solver
     virtual void backup_prevs();

     // pack and unpack the send/recv data for BD4 
    virtual void pack();   
    virtual void unpack(); 

    virtual size_t solnsize() { // the size of send/recv variables
        return dsize; 
    }
    
    virtual double* sendslns() {
        return sbuf; 
    }
     // unpack the received content into local variable(s) 
    virtual double* recvslns() {
        return rbuf;  
    }
     // write back the final solution (u_end) from correction item (Driver.pint_sum) 
    virtual void update_uend();
    
    virtual void rhs() {
        if(ndim==3)      rhs_ac_bd4_3d_(grid->nxyz, param.lamdaxyz, &nguard, b, soln, soln_1, G1, &dtk, &beta_);
        else if(ndim==2) rhs_ac_bd4_2d_(grid->nxyz, param.lamdaxyz, &nguard, b, soln, soln_1, G1, &dtk, &beta_);
        else if(ndim==1) rhs_ac_bd4_1d_(grid->nxyz, param.lamdaxyz, &nguard, b, soln, soln_1, G1, &dtk, &beta_);
    }

    virtual void stencil() {
        if(ndim==3)      stencil_ac_bd4_3d_(grid->nxyz, param.lamdaxyz, &nguard, A, soln, &dtk, &beta_);
        else if(ndim==2) stencil_ac_bd4_2d_(grid->nxyz, param.lamdaxyz, &nguard, A, soln, &dtk, &beta_);
        else if(ndim==1) stencil_ac_bd4_1d_(grid->nxyz, param.lamdaxyz, &nguard, A, soln, &dtk, &beta_);
    }

    virtual void rhs_g1() {
        if(ndim==3)      rhs_g1_ac_bd4_3d_(grid->nxyz, &nguard, soln_1, soln_2, soln_3, soln_4, G1, &dtk, &beta_);
        else if(ndim==2) rhs_g1_ac_bd4_2d_(grid->nxyz, &nguard, soln_1, soln_2, soln_3, soln_4, G1, &dtk, &beta_);
        else if(ndim==1) rhs_g1_ac_bd4_1d_(grid->nxyz, &nguard, soln_1, soln_2, soln_3, soln_4, G1, &dtk, &beta_);
    }

};
#endif
