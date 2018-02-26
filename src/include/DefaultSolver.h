#ifndef _PinT_DEFAULTSOLVER_H_
#define _PinT_DEFAULTSOLVER_H_ 1

#include "Solver.h"

/*
 *  The class is a template of time integrator for linear system,
 *  For simple finite differentiation cases,like Crank-Nicolson method,
 *  three-variable is enough, that is there are only the A (stencil), the b (RHS) and the x (unknown/solution). 
 *  The class can be directly used without any modification.   
 *  See 'heat' for example.
 */

class DefaultSolver : public Solver {

protected:
    double *soln;   // the current solution, pointer to the grid->u_f/u_c

    virtual void stencil() = 0; 
    
    virtual void rhs() = 0;
public:
    DefaultSolver(PinT *c, Grid *g):Solver(c,g) {
        hypre = getLS(conf, grid);   
    }; 

    DefaultSolver(PinT *c, Grid *g, bool isFS) : Solver(c,g,isFS) {
        hypre = getLS(conf, grid);   
    }

    virtual unsigned long evolve();         // evolve over a time slice 
};

#endif
