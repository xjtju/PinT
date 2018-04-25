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
    if( (isFine) || (conf->num_std != 4)){
        Driver::Abort("BD4 is used only as coarse solver, and parameter num_std must be 4!\n");
    }
    dsize = grid->size * conf->num_std;  // 4th order backward differentiation
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
    sbuf   = alloc_mem(dsize);
    rbuf   = alloc_mem(dsize);
    prevslns = alloc_mem(dsize);
    slns     = alloc_mem(dsize);
}

void BD4Solver::pack()   { // pack is easy
    //soln_3/2/1/0 -> sbuf
    bd4_pack_(sbuf, soln_3, soln_2, soln_1, soln, &size);
}

void BD4Solver::unpack() {
    // rbuf -> soln_1/2/3/4
    bd4_unpack_(rbuf, soln_4, soln_3, soln_2, soln_1, &size);
    // use the newest value as the initial guess for the staring point of the time slice    
    blas_cp_(grid->u_start, soln_1, &size);
}

double* BD4Solver::curr_solns(){
     bd4_pack_(slns, soln_3, soln_2, soln_1, soln, &size);
     return slns; 
}

void BD4Solver::backup_prevs() 
{
    blas_cp_(prevslns, this->curr_solns(), &dsize);  
    blas_cp_(grid->u_cprev, grid->u_c, &size);
}

void BD4Solver::update_uend() {
    bd4_update_uend_(grid->u_end, this->sendslns(), &size);
}

// initialization for the first time slice is different for the other slices
// 1     : soln1/2/3/4=soln
// others: unpack
void BD4Solver::init_holder() {
    if(conf->mytid == 0) {
        blas_cp_(soln_4, soln, &size); 
        blas_cp_(soln_3, soln, &size); 
        blas_cp_(soln_2, soln, &size); 
        blas_cp_(soln_1, soln, &size);
    }
//    printf("mytid:%d, bd4 solns: %e, %e, %e, %e \n",conf->mytid, soln_1[size/4], soln_2[size/4], soln_3[size/4], soln_4[size/4]);
}

void BD4Solver::update_holder(){
    blas_cp_(soln_4, soln_3, &size); 
    blas_cp_(soln_3, soln_2, &size); 
    blas_cp_(soln_2, soln_1, &size); 
    blas_cp_(soln_1, soln,   &size);
}
