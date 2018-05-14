#include "NewtonSolver.h"
#include "Driver.h"

// overwrite the default evolve for New-Raphson method
unsigned long NewtonSolver::evolve() {
    // step0: set initial value
    soln = getSoln();        // pointer to the start point  
    grid->guardcell(soln);   // make sure the guardcell is synchonized  

    init_holder();        // init holders from latest solution at the start of the current time slice
    blas_clear_(unk, &size);

    unsigned counter = 0;
    int iter = 0;
    for(int i=0; i<steps; i++){
        // step1 : set boundary condition
        grid->bc(soln); 
        // step2 : call the solver
        iter = newton_raphson();
        counter = counter + iter;
    }
    // step3: return latest solution to PinT framework 
    // nothing need to do 
    return counter;
}

// Newton-Raphson's iteration loop
int NewtonSolver::newton_raphson() {
    
    bool ifg = false; // converge flag
    double err = 0;   // eps check 
    double itr_eps;

    // step0 : set F_{n-1} and calcaluate RHS G1 
    update_holder();
    rhs_g1();

    int counter = 0;
    int iter = 0; 

    for(int i=0; i<newton_itmax; i++) {
        // step1 : set initial guess value 
        // the initial guess of unk must be consistent with the value of soln 
        blas_clear_(unk, &size);
     
        // step2 : calcaluate RHS: b 
        rhs();

        // step3 : set the stencil struct matrix 
        stencil();

        // step4 : call the linear solver
        if(ls_eps_dynamic) { // decrease the eps of linear solver  gradually
            itr_eps = 1.0/pow(ls_eps_factor, i);
            if(ls_eps < itr_eps)
                hypre->set_eps(itr_eps);
        }
        iter = hypre->solve(unk, b, A);  
        counter = counter + iter; 

        // step5: update solution 
        update();
        // when solution is changed, synchonization is necessary
        grid->guardcell(soln);
        // the boundary condition should be applied again 
        // NOTE : 
        //   sometimes the bc has a big impact on the iteration number of linear solver 
        grid->bc(soln);       

        // step5: check converge in all grids   
        chk_eps(&err);
        // NOTE : in parallel version, BREAK must be carefully used.
        // make sure all the MPI process will be jump out of the loop at the same point
        if(err <= newton_eps) { 
           ifg = true;  
           break;
        }
    }
    if(!ifg)
    {
        double s1 = 1.0;
        double u1 = 1.0; 
        // for debug
        blas_dot_3d_(grid->nxyz, &nguard, unk,  unk,  &u1 );
        blas_dot_3d_(grid->nxyz, &nguard, soln, soln, &s1 );
        Driver::Abort("Newton Raphson loop does not converge, eps=%e, |unk|=%e, |sln|=%e\n" , err, u1, s1);
    }
    return counter;
}

//NOTE : for space division, converge check must be performed in the whole geographical space 
void NewtonSolver::chk_eps(double *err) {
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
