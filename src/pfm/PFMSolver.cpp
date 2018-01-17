#include "PFMSolver.h"

/**
 * at current, PFM cannot be able to supported by Driver. 
 * NEED IMPROVED!
 */

PFMSolver::PFMSolver(PinT *c, Grid *g) : PBiCGStab(c,g){
    setup();
}

PFMSolver::PFMSolver(PinT *c, Grid *g, bool isFS) : PBiCGStab(c,g,isFS){
    setup();
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
void PFMSolver::setup(){
    this->eps = 1.0e-6;
    this->itmax = 10;
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

    F_ = alloc_mem(this->size);
    F  = alloc_mem(this->size);
    G1 = alloc_mem(this->size);

    for(int i=nguard; i<nx+nguard; i++){
        double x = grid->getX(i); 
        F[i] = 1.0 + tanh( (x-0.5)/xi);    // initial value 
        F[i] = 1.0 - 0.5*F[i];
        grid->set_val4all(i,F[i]);
    }
    blas_cp_(Fm, F, &size);  
}

double* PFMSolver::fetch() {
    return unk;
}
void PFMSolver::update() {
}

void PFMSolver::prepare() {
}

// overwrite the default evolve for New-Raphson method
void PFMSolver::evolve() {
     
     // step0: set initial value
    blas_cp_(F, grid->u_start, &size); 

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
    blas_cp_(unk, F, &size);  
    
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
        // step2 : set RHS cg_b1d 
        
        // ste3 : Call linear solver
        solve();  
        
        
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

// b = -F^{k-1} 
void PFMSolver::cg_b1d(double *x){
    double g2; 
    for(int i=nguard; i<nx+nguard; i++){
       g2 = lamda_x * ( F[i-1] - 2*F[i] + F[i+1] )
           - dtk * F[i] * (F[i]-1.0) * (F[i] - beta_);
       b[i] = - ( F[i] - F_[i] - theta*g2 - (1-theta)*G1[i] ) ;
    }
}

// matrix * vector,  the stencil matrix 
// v = Ay
void PFMSolver::cg_Xv1d(double* v, double *x) {
    for(int i=nguard; i<nx+nguard; i++){
       A0 = 1 + 2*theta*lamda_x + theta*dtk * ( 
             (F[i]-1.0) * (F[i] - beta_) 
           + (F[i]    ) * (F[i] - beta_)
           + (F[i]    ) * (F[i] - 1.0  )  );

       v[i] = A1_*x[i-1] + A0*x[i] + A1*x[i+1];
    }
}

//calcaluate the residual r = b - Ax
void PFMSolver::cg_rk1d(double *r, double *x, double *b){
    double val = 0.0;
    for(int i=nguard; i<nx+nguard; i++){
       A0 = 1 + 2*theta*lamda_x + theta*dtk * ( 
             (F[i]-1.0) * (F[i] - beta_) 
           + (F[i]    ) * (F[i] - beta_)
           + (F[i]    ) * (F[i] - 1.0  )  );
       val = A1_*x[i-1] + A0*x[i] + A1*x[i+1];
       r[i] = b[i] - val;
    }
}


//
// 2D, used Fortran
//

void PFMSolver::cg_rk2d(double *r, double *x, double *b){
}

// matrix * vector,  the stencil matrix 
// v = Ay
void PFMSolver::cg_Xv2d(double* v, double *y) {
}

void PFMSolver::cg_b2d(double *x){
}

//
// 3D, used Fortran
//
void PFMSolver::cg_rk3d(double *r, double *x, double *b){
}

// matrix * vector,  the stencil matrix 
// v = Ay
void PFMSolver::cg_Xv3d(double* v, double *y) {
}

void PFMSolver::cg_b3d(double *x){
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

