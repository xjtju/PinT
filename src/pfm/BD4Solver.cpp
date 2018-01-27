#include "BD4Solver.h"

BD4Solver::BD4Solver(PinT *c, Grid *g) :NewtonSolver(c, g) {
    param.init(c,g);
    setup(); 
} 

BD4Solver::BD4Solver(PinT *c, Grid *g, bool isFS) : NewtonSolver(c,g,isFS) {
    param.init(c, g, isFS);
    setup(); 

    if(0 == grid->myid)
        if(isFS) param.printLamda("Fine   Solver"); 
        else     param.printLamda("Coarse Solver"); 
} 

void BD4Solver::setup() {
    dtk = param.dtk;
    beta_ = param.beta_;
    create_holder();
}

/**
 * for code clear, the holders for structure perservation are created by an indepent function 
 */
void BD4Solver::create_holder() {
    //soln_1 has been created by superclass
    soln_2 = alloc_mem(size);
    soln_3 = alloc_mem(size);
    soln_4 = alloc_mem(size);
}

void BD4Solver::init_holder() {
    blas_cp_(soln_1, soln, &size); 
    blas_cp_(soln_2, soln, &size); 
    blas_cp_(soln_3, soln, &size); 
    blas_cp_(soln_4, soln, &size); 
}

void BD4Solver::update_holder(){
    blas_cp_(soln_4, soln_3, &size); 
    blas_cp_(soln_3, soln_2, &size); 
    blas_cp_(soln_2, soln_1, &size); 
    blas_cp_(soln_1, soln,   &size);
}
