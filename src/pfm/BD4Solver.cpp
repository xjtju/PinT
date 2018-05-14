#include "BD4Solver.h"

BD4Solver::BD4Solver(PinT *c, Grid *g) :NewtonSolver(c, g) {
    param.init(c,g);
    setup(); 
} 

BD4Solver::BD4Solver(PinT *c, Grid *g, bool isFS) : NewtonSolver(c,g,isFS) {
    param.init(c, g, isFS);
    setup(); 

    if(0 == grid->myid)
        if(isFS) param.printLamda("BD Fine Solver"); 
        else     param.printLamda("BD Coar Solver"); 
} 

void BD4Solver::setup() {
    if( (isFine) || (conf->num_std != 4)){
        Driver::Abort("BD4 is used only as coarse solver, and parameter 'num_std' must be 4!\n");
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
    bd4_pack_(sbuf, soln_3, soln_2, soln_1, getSoln(), &size);
}

/**
 * NOTE : ??? 
 *   the init of BD4 at the start of time slice has a big impact on convergence rate
 *     1: simplified BD4, soln_4 = soln_3 = soln_2 = soln_1 , faster convergence rate
 *     2. a variant BD4 , soln_4 = solon_4, soln_3 = soln_3, soln_2 = soln_1, 
 *     3: standard BD4,   soln_4 = soln_4, soln_3=soln_3, soln_2=soln_2, soln_1=soln_1 
 *
 *   According the 1D test case (k=1600, beta=-0.128, T=0.1, rfc_=1000) 
 *     NO.2 has the best convergence rate, 
 *     NO.1 is always faster one iteration than NO.3  
 *   
 */
void BD4Solver::unpack() {
    // rbuf -> soln_1/2/3/4
    bd4_unpack_(rbuf, soln_4, soln_3, soln_2, soln_1, &size);
    //blas_cp_(soln_4, soln_1, &size);
    //blas_cp_(soln_3, soln_1, &size);
    //blas_cp_(soln_2, soln_1, &size);
    // use the newest value as the initial guess for the staring point of the time slice    
    blas_cp_(grid->u_start, soln_1, &size);

    //double d4, d3, d2, d1;
    //blas_vdist_1d_(grid->nxyz, &grid->nguard, soln_4,   soln_1, &d4); 
    //blas_vdist_1d_(grid->nxyz, &grid->nguard, soln_3,   soln_1, &d3); 
    //blas_vdist_1d_(grid->nxyz, &grid->nguard, soln_2,   soln_1, &d2); 
    //blas_vdist_1d_(grid->nxyz, &grid->nguard, soln_1,   soln,   &d1); 
    //printf("mytid=%d, d1-4=%e, d1-3=%e, d1-2=%e, d0-1=%e\n", conf->mytid, d4, d3, d2, d1);
}

double* BD4Solver::curr_solns(){
    bd4_pack_(slns, soln_3, soln_2, soln_1, getSoln(), &size);
    return slns; 
}

void BD4Solver::backup_prevs() 
{
    blas_cp_(prevslns, this->curr_solns(), &dsize);  
    blas_cp_(grid->u_cprev, getSoln(), &size);
}

void BD4Solver::update_uend() {
    bd4_update_uend_(grid->u_end, this->sendslns(), &size);
}

// initialization for the first time slice is different for the other slices
// 1     : soln_1/2/3/4 <- soln or soln_1 <- soln
// others: unpack, nothing need to do here
void BD4Solver::init_holder() {
    if(conf->mytid == 0) {
        //blas_cp_(soln_4, getSoln(), &size); 
        //blas_cp_(soln_3, getSoln(), &size); 
        //blas_cp_(soln_2, getSoln(), &size); 
        blas_clear_(soln_4, &size);
        blas_clear_(soln_3, &size);
        blas_clear_(soln_2, &size);
        blas_cp_(soln_1, getSoln(), &size);
    }
}

// update four solns in order 
void BD4Solver::update_holder(){
    blas_cp_(soln_4, soln_3, &size); 
    blas_cp_(soln_3, soln_2, &size); 
    blas_cp_(soln_2, soln_1, &size); 
    blas_cp_(soln_1, getSoln(),   &size);

}
