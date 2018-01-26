#include "EulerSolver.h"

/**
 * For simple test only, not well tuned.
 * Euler method is not often used, except as a simple example of numerical integration.
 */ 
EulerSolver::EulerSolver(PinT *c, Grid *g) : PFMSolver(c,g){
}

EulerSolver::EulerSolver(PinT *c, Grid *g, bool isFS) : PFMSolver(c,g,isFS){
}
    
/**
 * Forward Euler, 1st order
 */
void EulerSolver::evolve() {

    // step0: set initial value
    soln = getSoln();     // pointer to the start point  

    for(int i=0; i<steps; i++){
        // step1 : set boundary condition
        grid->bc(soln); 
        // step2 : call the solver
        euler(); 
    }
    // step3: return latest solution to PinT framework 
    // nothing need to do 
}

void EulerSolver::euler() {

    rhs();    // calcaluate RHS 
    update(); // Xn+1 = Xn + RHS (delta) 

}
