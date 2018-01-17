#include "PFMSolver.h"

/**
 * at current, PFM cannot be able to supported by Driver. 
 * NEED IMPROVED!
 */

PFMSolver::PFMSolver(PinT *c, Grid *g) : Solver(c,g){
    setup(c,g);
}

PFMSolver::PFMSolver(PinT *c, Grid *g, bool isFS) : Solver(c,g,isFS){
    setup(c,g);

    if(isFS) {
        lamda_x = d*conf->f_dt/(grid->dx*grid->dx);   
        dtk = conf->f_dt*k;
    }else {
        lamda_x = d*conf->c_dt/(grid->dx*grid->dx);   
        dtk = conf->c_dt*k; 
    }

    A1_ = -theta*lamda_x;
    A1  = -theta*lamda_x;

    if(0 == grid->myid)
    printf("Solver (Fine=%d) : lamda_x=%f, lamday=%f, lamdaz=%f \n", isFS, lamda_x, lamda_x, lamda_x);
}

// set diffuse coefficient and tune the default parameter, problem specific
void PFMSolver::setup(PinT *c, Grid *g){
    ls_eps = 1.0e-6;
    ls_itmax = 10;
    //this->steps = 6000;
    conf->init_module(this, pfm_inih);
    if(grid->myid == 0) {
        printf("PFM init parameter : \n");  
        printf("        D    = %f \n", d);
        printf("        k    = %f \n", k);
        printf("        beta = %f \n", beta);
        printf("  newton_iter=%d \n\n", newton_itmax);
    }

    double xi = sqrt(2.0*d/k);
    beta_ = 0.5 - beta;

    unk = alloc_mem(this->size);
    F_  = alloc_mem(this->size);
    G1  = alloc_mem(this->size);
    b   = alloc_mem(this->size);
    if(ndim==1)
        bcp = alloc_mem(3*this->size); 

    double val;
    for(int i=nguard; i<nx+nguard; i++){
        double x = grid->getX(i); 
        val = 1.0 + tanh( (x-0.5)/xi);    // initial value 
        val = 1.0 - 0.5*val;
        grid->set_val4all(i,val);
    }

    hypre = new PBiCGStab(c, g);
}

// overwrite the default evolve for New-Raphson method
void PFMSolver::evolve() {
     
    // step0: set initial value
    F = getSoln();     // pointer to the start point  
    blas_clear_(unk, &size);

    for(int i=0; i<steps; i++){
    //   grid->bc(F_);
    // step1 : set boundary condition
        F[0] = 2.0 - F[1];
        F[nx+nguard] = -F[nx];

        // step2 : call the solver
        newton_raphson();
    }

    // step3: return solution 
    // nothing need to do 
}

void PFMSolver::newton_raphson() {
    
    blas_cp_(F_, F, &size); 
    for(int i=nguard; i<nx+nguard; i++){
        G1[i] = lamda_x * ( F_[i-1] - 2*F_[i] + F_[i+1] )
           - dtk * F_[i] * (F_[i]-1.0) * (F_[i] - beta_);
    }
    bool ifg = false;

    for(int i=0; i<newton_itmax; i++) {
        // step1 : set initial guess value
        blas_clear_(unk, &size);
        // step2 : set RHS 
        this->rhs();

        this->stencil();

        // ste3 : Call linear solver
        hypre->solve(unk, b, bcp);  
        
        // ste4: update solution 
        for(int j=nguard; j<nx+nguard; j++){
            F[j] = F[j] + unk[j];
        }

        // step5: check converge 
        double err =0.0;
        blas_dot_1d_(grid->nxyz, &nguard, unk, unk, &err ); 
        err = sqrt(err);
        if(err <= 1.0e-6) {
           ifg = true;
            break;
        }
    }
    if(!ifg) Driver::Abort("Newton Raphson loop does not converge");
}

//
// 1D, not used Fortran

void PFMSolver::rhs(){
    double g2; 
    for(int i=nguard; i<nx+nguard; i++){
       g2 = lamda_x * ( F[i-1] - 2*F[i] + F[i+1] )
           - dtk * F[i] * (F[i]-1.0) * (F[i] - beta_);
       b[i] = - ( F[i] - F_[i] - theta*g2 - (1-theta)*G1[i] ) ;
    }
}

/*
void PFMSolver::stencil(){
    int ind;
    for(int i=nguard; i<nx+nguard; i++){
       A0 = 1 + 2*theta*lamda_x + theta*dtk * ( 
             (F[i]-1.0) * (F[i] - beta_) 
           + (F[i]    ) * (F[i] - beta_)
           + (F[i]    ) * (F[i] - 1.0  )  );
    }
}
*/
// parse ini parameters of PFM
int pfm_inih(void* obj, const char* section, const char* name, const char* value) {

    PFMSolver* pfm = (PFMSolver*)obj;
    if (MATCH("pfm","ac_kval"))        { pfm->k = atof(value); }
    if (MATCH("pfm","ac_dval"))        { pfm->d = atof(value); }
    if (MATCH("pfm","ac_beta"))        { pfm->beta = atof(value); }
    if (MATCH("pfm", "newton_itmax"))  { pfm->newton_itmax = atoi(value); }

    return 0;
}

