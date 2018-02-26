// PBiCBSTAB 
#include "PBiCGStab.h"
#include "Driver.h"

/**
 * For the Heat equation with a constant diffuse factor or other simple examples, 
 * the coefficients of each stencil for different grid point may be same at the entire grid during the whole execution, 
 * so it is not necessary to use a LARGE matrix to hold these coefficients. 
 *
 * BUT in practice, for real world simulations, 
 * the coefficients associated with each stencil entry will typically vary from gridpoint to gridpoint,
 * so the caller must provide both the RHS(b) and stencil struct matrix(A) 
 *
 */

void PBiCGStab::init() {
    //b   = alloc_mem(size);

    r0_ = alloc_mem(size);
    r   = alloc_mem(size);
    v   = alloc_mem(size);
    t   = alloc_mem(size);

    p   = alloc_mem(size);
    p_  = alloc_mem(size);
    s   = alloc_mem(size);
    s_  = alloc_mem(size);
}

int PBiCGStab::solve(double *x, double *b, double *A){

    double rho0, rho, alpha, omega, beta;
    double res = 1000.0;
   
    rho0 = alpha = omega = 1.0;     
      
    blas_clear_(v,  &size);
    blas_clear_(r,  &size);
    blas_clear_(p,  &size);
    blas_clear_(p_, &size);

    grid->guardcell(x);    
    //cg_b(x);  // Ax=b, prepare b

    cg_rk(r, x, b, A); //init residual r = b - Ax 

    // choose an aribtrary vector r0_ such that (r0_, r) /=0, e.g. r0_ = r
    blas_cp_(r0_, r, &size);

    double tmp, tmp1, tmp2; 
    tmp = tmp1 = tmp2 = 0.0;  

    double b_nrm2 = cg_vdot(b, b);
    grid->sp_allreduce(&b_nrm2);
    int i;
    for(i=1; i<=itmax; i++){
        rho = cg_vdot(r0_, r);
        grid->sp_allreduce(&rho);
        //if( i==1 ){
        //   blas_cp_(p, r, &size);
        //}else {
        beta = (rho/rho0) * (alpha/omega);
        cg_direct(p, r, v, beta, omega); 
       //}

        grid->guardcell(p);
        preconditioner(p_, p, A, isPrecond); //solve Mp_=p, in some algorithm description, p_ is also denoted by y
        if(isPrecond) grid->guardcell(p_);
        cg_ax(v,p_,A);        // v=Ap_    

        tmp = cg_vdot(r0_, v);
        grid->sp_allreduce(&tmp);
        alpha = rho / tmp; 
        // h = x + alpha*p_, is ||h|| is sufficiently small, can do...  
        
        cg_avpy(s, alpha, r, v);  // s = r -alpha*v;
        grid->guardcell(s);
        // solve Ms_=s , in some algorithm description, s_ is also denoted by z
        preconditioner(s_, s, A, isPrecond); 
        
        if(isPrecond) grid->guardcell(s_);
        cg_ax(t,s_,A); // t=Az

        // in preconditioner t = 1/M_1*t and s = 1/M_1*s
        tmp1 = cg_vdot(t, s);
        grid->sp_allreduce(&tmp1);
        tmp2 = cg_vdot(t, t); // t dot t 
        grid->sp_allreduce(&tmp2);
        omega = tmp1 / tmp2;   
        
        cg_xi(x, alpha, p_, omega, s_); // xi = xi_1 + alpha*y + omega*z;   

        cg_avpy(r, omega, s, t);  // r = s - omega*t
         
        // ||r|| 
        tmp = cg_vdot(r, r);
        grid->sp_allreduce(&tmp);   
        res = sqrt(tmp/b_nrm2);

        if(res < eps) {
            break;
        }
        rho0 = rho;
        // WARN: the bc maybe introduce numerical error when grid default bc function is not used for some problems
        // In this condition, problem-specific sub grid must be implemented
        // grid->bc(x); // necessary ?
    }
    if(i>=itmax) Driver::Abort("PBiCG is not converged: Iter: %d, rho : %e, beta : %e, alpha : %e, omega : %e, res : %e \n",i,  rho, beta, alpha, omega, res);

    return i;
    //grid->guardcell(x); // necessary ?
}

//
// vector operations can be accelerated by SIMD. 
// NOTE: the following 1d functions are for test only, in practice, their fortran counterparts are used instead 
//

//p = r + beta * ( p - omg * v )
void PBiCGStab::cg_direct1d(double* p, double* r, double* v, double beta, double omega){
    for(int i=nguard; i<nx+nguard; i++) {
        p[i] = r[i] + beta*(p[i] - omega*v[i]);
    }
}

// x = x + ay + bz  
void PBiCGStab::cg_xi1d(double* x, double* y, double* z, double alpha, double omega){
    for(int i=nguard; i<nx+nguard; i++) {
        x[i] = x[i] + alpha*y[i] + omega*z[i];
    }
}

//
// 2D, 
// NOTE: the following 2d functions are for test only, in practice, their fortran counterparts are used instead 
//

void PBiCGStab::cg_direct2d(double* p, double* r, double* v, double beta, double omega){
    long idx;
    for(int j=nguard; j<ny+nguard; j++)
        for(int i=nguard; i<nx+nguard; i++) {
            idx = grid->getOuterIdx(i,j,0);
            p[idx] = r[idx] + beta*(p[idx] - omega*v[idx]);
    }
}

void PBiCGStab::cg_xi2d(double* x, double* y, double* z, double alpha, double omega){
    long idx;
    for(int j=nguard; j<ny+nguard; j++) 
        for(int i=nguard; i<nx+nguard; i++) {
            idx = grid->getOuterIdx(i,j,0); 
            x[idx] = x[idx] + alpha*y[idx] + omega*z[idx];
    }
}


/**
 * preconditioner, in our original version ,the below idea was adapted. 
 * solve Mp_=p
 * M = M1*M2 ~ A, there is a simple method for decomposing the matrix A when A has some regular pattern, 
 * for example in HEAT diffusion:
 * if we set  M1=I, M2=A, where I is Identity matrix.  
 * then the preconditioner can be set to a light-weight but less precise solver than BiCG, 
 * and solving the same linear system: Ax=b, moreover the following calculations of BiCG will become : 
    ==  reverse(M1)*t=I*t=t and reverse(M1)*s=I*s=s ==  
 * so the matrix reverse and matrix multiplying vector calculations can also be skipped.   
 *
 * BUT: the above approach did not worked, that is NO EFFECT, moreover, when iterations < 10, even become worse. 
 * so we have to do matrix decomposition, and now there is no effective preconditioner in the current version.
 */
void PBiCGStab::preconditioner(double* p_, double* p, double *A, bool isPrecond){
    if(!isPrecond)   
        blas_cp_(p_, p, &size);
    else {  
        sor->solve(p_, p, A);
    }
}

