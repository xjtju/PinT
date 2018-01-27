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
    
    theta = param.theta;  
    beta_ = param.beta_;       
    dtk   = param.dtk;         

    // overwrite the default value from .ini file
    this->newton_itmax = param.newton_itmax;

    //this->steps = 1; // only for debug
    //conf->kpar_limit = 0; //only for debug
}
