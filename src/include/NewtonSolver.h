#ifndef _PinT_NEWTONSOLVER_H_
#define _PinT_NEWTONSOLVER_H_ 1

#include "Solver.h"

/**
 *
 *  The class is a template of Newton-Raphson method for nonlinear system 
 *  
 *  For example, the typical phase field model, Allen-Cahn equation  
 *  Compared to the HEAT, it is introduced extra three variables (G1, soln_1, unk) for data structure preservation 
 *  during the calculation due to its nonlinear feature. 
 *
 *  Proper differentiation formulas can be easily implemented by overwriting all the virtual functions.
 *  See the CNSolver (Crank-Nicolson) and BD4Solver (4th order backward euler) for example in PFM.
 *
 */

class NewtonSolver : public Solver {

private:
    bool ls_eps_dynamic;
    double ls_eps;
    double ls_eps_factor;

    void setup(){
        ls_eps_dynamic = conf->ls_eps_dynamic;
        ls_eps_factor = conf->ls_eps_factor;
        ls_eps = conf->ls_eps;

        unk   = alloc_mem(this->size);
        soln_1= alloc_mem(this->size);
        G1    = alloc_mem(this->size);
        hypre = getLS(conf, grid);   // any linear solver implementing LS is OK
    }

public:
    NewtonSolver(PinT *c, Grid *g):Solver(c,g) {
        setup(); 
    }; 

    NewtonSolver(PinT *c, Grid *g, bool isFS) : Solver(c,g,isFS) {
        setup(); 
    }

    virtual ~NewtonSolver() {

        free_mem(soln_1);
        free_mem(G1);
        free_mem(unk);

        if(grid->myid==0 && conf->verbose)
        printf("INFO: the memory allocate by NewtonSolver has been released.\n");
    }

    virtual unsigned long evolve();         // evolve over a time slice 

protected:
    int newton_itmax = 50;
    double newton_eps= 1.0e-6;

    double *soln;   // the current solution, pointer to the grid->u_f/u_c
    // structure preservation 
    double *soln_1;  // the holder of -F^{k-1} in Newton's method when applying to nonlinear systems of equations 

    /** 
     * the U^{n-1} part of RHS,
     *    CN : F^{n} = U^{n}-U^{n-1} - theta*G2 - (1-theta)*G1 
     *   BD4 : F^{n} = 25*U^{n} + G1 - 12*G2
     * Because there are iterations in Newton method, the G1 is calcaluated just once before entering iterations, 
     * but G2 is necessary to be re-calcaluated in each iteration, so G1 is reserved for reuse, and G2 is temporary.     
     */
    double *G1;      
    double *unk;     // the unknown 'x' of Ax=b, that is (U_{n+1} - U_{n}) in Newton's method 

    int newton_raphson(); // the template of New-Raphson method iteration 
    
    /*
     * NOTE : for space division, converge check must be performed in the whole geographical space 
     */
    void chk_eps(double *err); 

    // init previous solution holder for backward method, the structure preservation
    // at the starting, all of them are identical to the initial condition of the problem 
    virtual void init_holder() {
        blas_cp_(soln_1, soln, &size); 
    }

    // update backward solution holders for the next step 
    virtual void update_holder(){
        blas_cp_(soln_1, soln, &size); 
    }
    
    virtual void stencil() = 0; 
    
    virtual void rhs() = 0;   // calcaluate the whole RHS
    
    virtual void rhs_g1() = 0; // calcaluate G1 part of RHS 

    void update() {
        if(ndim==3)      update_soln_3d_(grid->nxyz, &nguard, soln, unk);
        else if(ndim==2) update_soln_2d_(grid->nxyz, &nguard, soln, unk);
        else if(ndim==1) update_soln_1d_(grid->nxyz, &nguard, soln, unk);
    }

};

#endif
