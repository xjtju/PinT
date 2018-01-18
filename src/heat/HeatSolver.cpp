#include "HeatSolver.h"

HeatSolver::HeatSolver(PinT *c, Grid *g) : Solver(c,g){
    setup();
}

HeatSolver::HeatSolver(PinT *c, Grid *g, bool isFS) : Solver(c,g,isFS){
    setup();
}

// set diffuse coefficient and tune the default parameter, problem specific
void HeatSolver::setup(){
    if(ndim==1) k = 0.061644; 
    if(ndim==2) k = 0.061644;
    if(ndim==3) k = 0.061644;

    hypre = new PBiCGStab(conf, grid); // choose a linear solver
}

// evolve along one time slice for Crank-Nicolson method  
void HeatSolver::evolve() {
     
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
        hypre->solve(soln, b, bcp);

        // step5: update solution, 
        // for heat diffusion, soln has already updated by linear solver, noting need to be done 
    }
}
