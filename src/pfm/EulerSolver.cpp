#include "EulerSolver.h"

/**
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
        if(isFS) param.printLamda("EU Fine Solver"); 
        else fprintf(stderr, "WARN : EulerSolver is used as coarse solver, unstable !\n"); 
}
    
void EulerSolver::setup(){
    beta_ = param.beta_;       
    dtk   = param.dtk;
    dsize = size*conf->num_std;
    rfc_ = conf->rfc_;
    if(conf->num_std == 4){
        soln_1 = alloc_mem(size);
        soln_2 = alloc_mem(size);
        soln_3 = alloc_mem(size);
        slns   = alloc_mem(dsize);
    }
}
/**
 * Forward Euler, 1st order
 */
unsigned long EulerSolver::evolve() {

    // step0: set initial value
    soln = getSoln();     // pointer to the start point  
    grid->guardcell(soln); // make sure the guardcell is at synchronization state  

    unsigned long counter = 0;
    for(long i=0; i<steps; i++){
        // step1 : set boundary condition
        grid->bc(soln); 
        // step2 : call the solver
        euler(i); 
        counter++;
        grid->guardcell(soln); // do not forget synchronize guardcell 
    }
    // step3: return latest solution to PinT framework 
    // nothing need to do 
    return counter;
}

// preserve the solutions at the correponding time steps of the coarse solver
void EulerSolver::euler(long step) {
    if(conf->num_std == 4) {
        if     ( step == (steps-3*rfc_) )   
            blas_cp_(soln_3, soln, &size); 
        else if( step == (steps-2*rfc_) )  
            blas_cp_(soln_2, soln, &size); 
        else if(  step == (steps-1*rfc_) ) 
            blas_cp_(soln_1, soln, &size); 
    }

    rhs();    // calcaluate RHS 
    update(); // Xn+1 = Xn + RHS (delta)
}

double* EulerSolver::curr_solns(){
    if(conf->num_std == 4) {
        bd4_pack_(slns, soln_3, soln_2, soln_1, getSoln(), &size);
        //double d4, d3, d2, d1;
        //blas_vdist_1d_(grid->nxyz, &grid->nguard, soln_3,   soln, &d4); 
        //blas_vdist_1d_(grid->nxyz, &grid->nguard, soln_2,   soln, &d3); 
        //blas_vdist_1d_(grid->nxyz, &grid->nguard, soln_1,   soln, &d2); 
        //printf("mytid=%d/%d, d3-0=%e, d2-0=%e, d1-0=%e\n", conf->mytid, rfc_, d4, d3, d2 );
        return slns;
    }
    else return getSoln();  
}
