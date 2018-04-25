#include "CNSolver.h"

/**
 *  Time integrator of Allen-Cahn Equation (AC)
 */
CNSolver::CNSolver(PinT *c, Grid *g) : NewtonSolver(c,g){
    param.init(conf, grid);
    setup();
}

CNSolver::CNSolver(PinT *c, Grid *g, bool isFS) : NewtonSolver(c,g,isFS){
    param.init(conf, grid, isFS);
    setup();

    if(0 == grid->myid)
        if(isFS) param.printLamda("Fine   Solver"); 
        else     param.printLamda("Coarse Solver"); 
}

// set diffuse coefficient and tune the default parameter, problem specific
void CNSolver::setup(){
    
    if((!isFine) && (conf->num_std!=1) ){
        Driver::Abort("parameter num_std must be 1 when CNSolver is as Coarse Solver!\n");
    }
    theta = param.theta;  
    beta_ = param.beta_;       
    dtk   = param.dtk;         

    dsize = size*conf->num_std;
    if( isFine && (conf->num_std == 4) ) {
        soln_1 = alloc_mem(size);
        soln_2 = alloc_mem(size);
        soln_3 = alloc_mem(size);
        slns   = alloc_mem(dsize);
    }

    // overwrite the default value from .ini file
    this->newton_itmax = param.newton_itmax;

    //this->steps = 1; // only for debug
    //conf->kpar_limit = 0; //only for debug
}

double* CNSolver::curr_solns(){
    if(conf->num_std == 4) {
        bd4_pack_(slns, soln_3, soln_2, soln_1, soln, &size);
        return slns;
    }
    else return getSoln();  
}
void CNSolver::update_holder(){
    if( isFine && (conf->num_std ==4) ) {
        blas_cp_(soln_3, soln_2, &size); 
        blas_cp_(soln_2, soln_1, &size); 
    }
    blas_cp_(soln_1, soln,   &size);
}
