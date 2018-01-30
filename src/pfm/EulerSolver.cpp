#include "EulerSolver.h"

/**
 * For simple test only, not well tuned.
 * Euler method is not often used, except as a simple example of numerical integration.
 */ 
EulerSolver::EulerSolver(PinT *c, Grid *g) : Solver(c,g)
{
    EulerSolver(c, g, true);
}

EulerSolver::EulerSolver(PinT *c, Grid *g, bool isFS) : Solver(c,g,isFS){
    param.init(c, g, isFS);
    setup();
    if(0 == grid->myid)
        if(isFS) param.printLamda("Fine   Solver"); 
        else fprintf(stderr, "WARN : EulerSolver is used as coarse solver, unstable !"); 
}
    
void EulerSolver::setup(){
    beta_ = param.beta_;       
    dtk   = param.dtk;         
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

        grid->guardcell(soln); // do not forget synchonize guardcell 
    }
    // step3: return latest solution to PinT framework 
    // nothing need to do 
}

void EulerSolver::euler() {

    rhs();    // calcaluate RHS 
    update(); // Xn+1 = Xn + RHS (delta) 
}
