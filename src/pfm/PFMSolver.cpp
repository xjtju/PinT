#include "PFMSolver.h"

/**
 *  Time integrator of Allen-Cahn Equation (AC)
 */
PFMSolver::PFMSolver(PinT *c, Grid *g) : Solver(c,g){
    setup();
}

PFMSolver::PFMSolver(PinT *c, Grid *g, bool isFS) : Solver(c,g,isFS){
    setup();

    if(isFS) {
        lamda_x = d*conf->f_dt/(grid->dx*grid->dx);   
        if(ndim>=2) lamda_y = d * conf->f_dt / (g->dy*g->dy);
        if(ndim>=3) lamda_z = d * conf->f_dt / (g->dz*g->dz);
        dtk = conf->f_dt*k;
    }else {
        lamda_x = d*conf->c_dt/(grid->dx*grid->dx);   
        if(ndim>=2) lamda_y = d * conf->c_dt / (g->dy*g->dy);
        if(ndim>=3) lamda_z = d * conf->c_dt / (g->dz*g->dz);
        dtk = conf->c_dt*k; 
    }

    lamdaxyz[0] = lamda_x;
    lamdaxyz[1] = lamda_y;
    lamdaxyz[2] = lamda_z;

    if(0 == grid->myid)
    printf("Solver (Fine=%d) : lamda_x=%f, lamday=%f, lamdaz=%f, dtk=%f\n", isFS, lamda_x, lamda_y, lamda_z, dtk);
}

// set diffuse coefficient and tune the default parameter, problem specific
void PFMSolver::setup(){
    conf->init_module(this, pfm_inih);
    beta_ = 0.5 - beta;
    if(grid->myid == 0) {
        printf("PFM init parameter : \n");  
        printf("        D    = %f \n", d);
        printf("        k    = %f \n", k);
        printf("        beta = %f \n", beta);
        printf("  newton_iter= %d \n", newton_itmax);
        printf("        theta= %f \n",theta);

        printf("     => 0.0: Ex.Euler; 0.5: Crank-Nicolson; 1.0: Im.Euler.\n\n"); 
    }

    unk   = alloc_mem(this->size);
    soln_ = alloc_mem(this->size);
    G1    = alloc_mem(this->size);

    hypre = getLS(conf,grid); 
    
    //this->steps = 1; // only for debug
    //conf->kpar_limit = 0; //only for debug
}

// overwrite the default evolve for New-Raphson method
void PFMSolver::evolve() {
    // step0: set initial value
    soln = getSoln();     // pointer to the start point  
    blas_clear_(unk, &size);

    for(int i=0; i<steps; i++){
        // step1 : set boundary condition
        grid->bc(soln); 
        // step2 : call the solver
        newton_raphson();
    }
    // step3: return latest solution to PinT framework 
    // nothing need to do 
}

// Newton-Raphson's iteration loop
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
        // when solution is changed, synchonization is necessary
        grid->guardcell(soln); 

        // step5: check converge in all grids   
        chk_eps(&err);
        // NOTE : in parallel version, BREAK must be carefully used.
        // make sure all the MPI process will be jump out of the loop at the same point
        if(err <= 1.0e-6) { 
           ifg = true;  
           break;
        }
    }
    if(!ifg) Driver::Abort("Newton Raphson loop does not converge, eps=%e\n" , err);
}

//NOTE : for space division, converge check must be performed in the whole geographical space 
void PFMSolver::chk_eps(double *err) {
        double x_nrm2 = 1.0;
        if(ndim==3) { 
            blas_dot_3d_(grid->nxyz, &nguard, unk,  unk,  err     );
            blas_dot_3d_(grid->nxyz, &nguard, soln, soln, &x_nrm2 );
        }
        else if(ndim==2) {
            blas_dot_2d_(grid->nxyz, &nguard, unk,  unk,  err     ); 
            blas_dot_2d_(grid->nxyz, &nguard, soln, soln, &x_nrm2 );
        }
        else if(ndim==1) {
            blas_dot_1d_(grid->nxyz, &nguard, unk,  unk,  err     );
            blas_dot_1d_(grid->nxyz, &nguard, soln, soln, &x_nrm2 );
        }

        grid->sp_allreduce(err);
        grid->sp_allreduce(&x_nrm2);
        *err = sqrt(*err/x_nrm2);
}

// set the initial value
void PFMSolver::init() {
    if(this->ndim == 1) init1d();
    if(this->ndim == 2) init2d();
    if(this->ndim == 3) init3d();
}
// set the initial value
void PFMSolver::init1d() {

    double xi = sqrt(2.0*d/k);
    double val;
    double midx = conf->Xspan/2 ; 
    double xdist;
    for(int i=nguard; i<nx+nguard; i++){
        double x = grid->getX(i); 
        val = 1.0 + tanh( (x-0.5)/xi);    // initial value 
        val = 1.0 - 0.5*val;
        //xdist = x - midx; 
        //if(abs(xdist) < 0.2) 
        //   val = 0.8;
        //else val = 0;
        grid->set_val4all(i,val);
    }
}

void PFMSolver::init2d(){
    long ind = 0;     
    double xdist, ydist, unk;
    for(int j = nguard; j<ny+nguard; j++)
    for(int i = nguard; i<nx+nguard; i++){
       xdist = grid->getX(i) -  conf->Xspan/2 ; 
       ydist = grid->getY(j) -  conf->Yspan/2 ;
       
       ind = grid->getOuterIdx(i, j, 0);  
       if( abs(xdist)<=0.4 &&  abs(ydist)<=0.4 )
           unk = 1.0; 
       else unk = 0.0; 
       
       grid->set_val4all(ind, unk);
    }
}

void PFMSolver::init3d(){
    long ind = 0;     
    double xdist, ydist, zdist, unk;
    for(int k = nguard; k<nz+nguard; k++)
    for(int j = nguard; j<ny+nguard; j++)
    for(int i = nguard; i<nx+nguard; i++){
       xdist = grid->getX(i) -  conf->Xspan/2 ; 
       ydist = grid->getY(j) -  conf->Yspan/2 ;
       zdist = grid->getZ(k) -  conf->Zspan/2 ;
       ind = grid->getOuterIdx(i, j, k);  
       if( (abs(xdist)<=0.4) &&  (abs(ydist)<=0.4) && (abs(zdist)<=0.4) )
           unk = 1.0; 
       else unk = 0.0; 
       
       grid->set_val4all(ind, unk);
    }
}
// parse ini parameters of PFM
int pfm_inih(void* obj, const char* section, const char* name, const char* value) {

    PFMSolver* pfm = (PFMSolver*)obj;
    if (MATCH("pfm","ac_kval"))        { pfm->k = atof(value); }
    if (MATCH("pfm","ac_dval"))        { pfm->d = atof(value); }
    if (MATCH("pfm","ac_beta"))        { pfm->beta = atof(value); }
    if (MATCH("pfm", "newton_itmax"))  { pfm->newton_itmax = atoi(value); }
    if (MATCH("pfm", "theta"))         { pfm->theta = atof(value); }

    return 0;
}

