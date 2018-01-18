#include "PFMSolver.h"

/**
 * at current, PFM cannot be able to supported by Driver. 
 * NEED IMPROVED!
 */

PFMSolver::PFMSolver(PinT *c, Grid *g) : Solver(c,g){
    setup();
}

PFMSolver::PFMSolver(PinT *c, Grid *g, bool isFS) : Solver(c,g,isFS){
    setup();

    if(isFS) {
        lamda_x = d*conf->f_dt/(grid->dx*grid->dx);   
        dtk = conf->f_dt*k;
    }else {
        lamda_x = d*conf->c_dt/(grid->dx*grid->dx);   
        dtk = conf->c_dt*k; 
    }

    if(0 == grid->myid)
    printf("Solver (Fine=%d) : lamda_x=%f, lamday=%f, lamdaz=%f, dtk=%f\n", isFS, lamda_x, lamda_x, lamda_x, dtk);
}

// set diffuse coefficient and tune the default parameter, problem specific
void PFMSolver::setup(){
    ls_eps = 1.0e-6;
    ls_itmax = 10;
    //this->steps = 1;
    conf->init_module(this, pfm_inih);
    if(grid->myid == 0) {
        printf("PFM init parameter : \n");  
        printf("        D    = %f \n", d);
        printf("        k    = %f \n", k);
        printf("        beta = %f \n", beta);
        printf("  newton_iter=%d \n\n", newton_itmax);
    }

    beta_ = 0.5 - beta;

    unk   = alloc_mem(this->size);
    soln_ = alloc_mem(this->size);
    G1    = alloc_mem(this->size);

    hypre = new PBiCGStab(conf, grid);
}

// set the initial value
void PFMSolver::init() {

    double xi = sqrt(2.0*d/k);
    double val;
    for(int i=nguard; i<nx+nguard; i++){
        double x = grid->getX(i); 
        val = 1.0 + tanh( (x-0.5)/xi);    // initial value 
        val = 1.0 - 0.5*val;
        grid->set_val4all(i,val);
    }

}
// overwrite the default evolve for New-Raphson method
void PFMSolver::evolve() {
     
    // step0: set initial value
    soln = getSoln();     // pointer to the start point  
    blas_clear_(unk, &size);

    for(int i=0; i<steps; i++){
        // step1 : set boundary condition
        bc();
        // step2 : call the solver
        newton_raphson();
    }

    // step3: return solution 
    // nothing need to do 
}

void PFMSolver::newton_raphson() {
    
    bool ifg = false; // converge flag
    double err = 0;   // eps check 
    
    // step0 : set F_{n-1} and calcaluate RHS G1 
    blas_cp_(soln_, soln, &size); 
    rhs_g1();

    for(int i=0; i<newton_itmax; i++) {
        // step1 : set initial guess value
        blas_clear_(unk, &size);

        // step2 : calcaluate RHS: b 
        rhs();
        
        // step3 : set the stencil struct matrix 
        stencil();

        // step4 : call the linear solver
        hypre->solve(unk, b, bcp);  
        
        // step5: update solution 
        update();

        // step5: check converge 
        chk_eps(&err);
        if(err <= 1.0e-6) {
           ifg = true;  break;
        }
    }
    if(!ifg) Driver::Abort("Newton Raphson loop does not converge, eps=%e\n" , err);
}

// parse ini parameters of PFM
int pfm_inih(void* obj, const char* section, const char* name, const char* value) {

    PFMSolver* pfm = (PFMSolver*)obj;
    if (MATCH("pfm","ac_kval"))        { pfm->k = atof(value); }
    if (MATCH("pfm","ac_dval"))        { pfm->d = atof(value); }
    if (MATCH("pfm","ac_beta"))        { pfm->beta = atof(value); }
    if (MATCH("pfm", "newton_itmax"))  { pfm->newton_itmax = atoi(value); }

    return 0;
}

