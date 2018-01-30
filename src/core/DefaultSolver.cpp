#include "DefaultSolver.h"

// the template for time integration of simple linear system 
// evolve along one time slice for Crank-Nicolson method  
void DefaultSolver::evolve() {
     
    // step0: set initial value, the latest solution is used directly 
    soln = getSoln();     // pointer to the start point of this time slice  

    for(int i=0; i<steps; i++){
        // step1 : set boundary condition, default bc function provided by Grid is enough
        grid->bc(soln);
        
        // step2 : calcaluate RHS
        rhs();

        // step3 : set the stencil struct matrix 
        stencil();

        // step4 : call the linear solver, not necessary to use a new guess for heat diffusion 
        hypre->solve(soln, b, A);

        // step5: update solution, 
        // for heat diffusion, soln has already updated by linear solver, noting need to be done 
        grid->guardcell(soln); 
    }
}
