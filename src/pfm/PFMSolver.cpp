#include "PFMSolver.h"

PFMSolver::PFMSolver(PinT *c, Grid *g) : PBiCGStab(c,g){
    setup();
}

PFMSolver::PFMSolver(PinT *c, Grid *g, bool isFS) : PBiCGStab(c,g,isFS){
    setup();
}

// set diffuse coefficient and tune the default parameter, problem specific
void PFMSolver::setup(){
    this->eps = 1.0e-6;
    conf->init_module(this, pfm_inih);

    if(grid->myid == 0) {
        printf("PFM init parameter : \n");  
        printf("  D   = %f \n", d);
        printf("  k   = %f \n", k);
        printf("  beta= %f \n\n", beta);
    }
}

//
// 1D, not used Fortran
//

//calcaluate the residual r = b - Ax
void PFMSolver::cg_rk1d(double *r, double *x, double *b){
    for(int i=nguard; i<nx+nguard; i++){

    }
}

// matrix * vector,  the stencil matrix 
// v = Ay
void PFMSolver::cg_Xv1d(double* v, double *y) {
    for(int i=nguard; i<nx+nguard; i++){

    }
}

void PFMSolver::cg_b1d(double *x){
    for(int i=nguard; i<nx+nguard; i++){

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

int pfm_inih(void* obj, const char* section, const char* name, const char* value) {

    PFMSolver* pfm = (PFMSolver*)obj;
    if (MATCH("pfm","ac_kval"))        { pfm->k = atof(value); }
    if (MATCH("pfm","ac_dval"))        { pfm->d = atof(value); }
    if (MATCH("pfm","ac_beta"))        { pfm->beta = atof(value); }

    return 0;
}

