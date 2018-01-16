#include "PFMSolver.h"

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
    this->eps = 1.0e-5;
    this->itmax = 10;
    //this->steps = 100;
    conf->init_module(this, pfm_inih);
    if(grid->myid == 0) {
        printf("PFM init parameter : \n");  
        printf("  D   = %f \n", d);
        printf("  k   = %f \n", k);
        printf("  beta= %f \n\n", beta);
    }

    beta_ = 0.5 - beta;
    unk = alloc_mem(this->size);
    F_ = alloc_mem(this->size);
    F  = alloc_mem(this->size);
    G1 = alloc_mem(this->size);

    prepare();
}

double* PFMSolver::fetch() {
   blas_clear_(unk, &size);  // init 0 
   return unk; 
}

void PFMSolver::update() {
    for(int i=nguard; i<nx+nguard; i++){
        double val;
        F[i] = F_[i] + unk[i];
        //blas_vdist_1d_(grid->nxyz, &grid->nguard, F, F_, &val); 
        F_[i] = F[i]; // structure preservation 
        G1[i] = lamda_x * ( F_[i-1] - 2*F_[i] + F_[i+1] )
           - dtk * F_[i] * (F_[i]-1.0) * (F_[i] - beta_);
        grid->u_end[i] = F[i];
    }
    //printf("\n\n");
}

void PFMSolver::prepare() {
    printf("preparing structure preservation\n");
    for(int i=nguard; i<nx+nguard; i++){
        F_[i] = grid->u[i]; 
        F[i]  = grid->u[i];
        G1[i] = lamda_x * ( F_[i-1] - 2*F_[i] + F_[i+1] )
           - dtk * F_[i] * (F_[i]-1.0) * (F_[i] - beta_);
    }
    //update();
}

// overwrite the default evolve for New-Raphson method
void PFMSolver::evolve() {
    
}
//
// 1D, not used Fortran

// b = -F^{k-1} 
void PFMSolver::cg_b1d(double *x){
    double g2; 
    for(int i=nguard; i<nx+nguard; i++){
       g2 = lamda_x * ( F[i-1] - 2*F[i] + F[i+1] )
           - dtk * F[i] * (F[i]-1.0) * (F[i] - beta_);
       b[i] = - ( (F[i] - F_[i]) - theta*g2 + (1-theta)*G1[i] ) ;
    }
}

// matrix * vector,  the stencil matrix 
// v = Ay
void PFMSolver::cg_Xv1d(double* v, double *x) {
    for(int i=nguard; i<nx+nguard; i++){
       A0 = 1 + 2*theta*lamda_x + theta*dtk * ( 
             (F_[i]-1.0) * (F_[i] - beta_) 
           + (F_[i]    ) * (F_[i] - beta_)
           + (F_[i]    ) * (F_[i] - 1.0  )  );

       v[i] = A1_*x[i-1] + A0*x[i] + A1*x[i+1];
    }
}

//calcaluate the residual r = b - Ax
void PFMSolver::cg_rk1d(double *r, double *x, double *b){
    double val = 0.0;
    for(int i=nguard; i<nx+nguard; i++){
       A0 = 1 + 2*theta*lamda_x + theta*dtk * ( 
             (F_[i]-1.0) * (F_[i] - beta_) 
           + (F_[i]    ) * (F_[i] - beta_)
           + (F_[i]    ) * (F_[i] - 1.0  )  );
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

    return 0;
}

