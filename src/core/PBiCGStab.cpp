// PBiCBSTAB 
#include "PBiCGStab.h"

void PBiCGStab::solve(){

    double *x = fetch();  // get the grid variables

    double rho0, rho, alpha, omega, beta;
    double res = 1000.0;
   
    rho0 = alpha = omega = 1.0;     
      
    blas_clear_(v,  &size);
    blas_clear_(b,  &size);
    blas_clear_(r,  &size);
    blas_clear_(p,  &size);
    blas_clear_(p_, &size);

    grid->guardcell(x);    
    cg_b(x);

    cg_rk(r, x, b); //init residual r = b - Ax 
    //grid->output_var(r,false);
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
        //   blas_cp(p, r, size);
        //}else {
            beta = (rho/rho0) * (alpha/omega);
            cg_direct(p, r, v, beta, omega); 
       //}

        grid->guardcell(p);
        preconditioner(p_, p, isPrecond); //solve Mp_=p, in some algorithm description, p_ is also denoted by y

        grid->guardcell(p_);
        cg_Xv(v,p_);        // v=Ap_    

        tmp = cg_vdot(r0_, v);
        grid->sp_allreduce(&tmp);
        alpha = rho / tmp; 
        // h = x + alpha*p_, is ||h|| is sufficiently small, can do...  
        
        cg_avpy(s, alpha, r, v);  // s = r -alpha*v;
        grid->guardcell(s);
        // solve Ms_=s , in some algorithm description, s_ is also denoted by z
        preconditioner(s_, s, isPrecond); 
        
        grid->guardcell(s_);
        cg_Xv(t,s_); // t=Az

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
        grid->bc(x);
    }

    //printf("Iter: %d, rho : %f, beta : %f, alpha : %f, omega : %f, res : %e \n",i,  rho, beta, alpha, omega, res);
    grid->guardcell(x);
    update(); // update grid variables
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
            idx = grid->getOuterIdx(i,j);
            p[idx] = r[idx] + beta*(p[idx] - omega*v[idx]);
    }
}

void PBiCGStab::cg_xi2d(double* x, double* y, double* z, double alpha, double omega){
    long idx;
    for(int j=nguard; j<ny+nguard; j++) 
        for(int i=nguard; i<nx+nguard; i++) {
            idx = grid->getOuterIdx(i,j); 
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
void PBiCGStab::preconditioner(double* p_, double* p, bool isPrecond){
    if(!isPrecond)   
        blas_cp_(p_, p, &size);
    else {  
        sor2(p_, p, 10, false);
    }
}

/**
 * check convergency flag is not used. if used, the previous solution must be hold for comparing. 
 */
void PBiCGStab::sor2(double *p_, double *p, int lc_max, bool checkCnvg){
    int lc = 0;
    double res = 0.0;

    for(lc=1; lc<=lc_max; lc++) {
        for(int color=0; color<2; color++) {
            grid->guardcell(p_);
            sor2_core(p_, p, &color);
         }
        if(checkCnvg) {
            if( res < eps) break;
        }
    }
}

void PBiCGStab::update() {
    //in current implementation, the solver directly modify the grid variables, 
    //so it is not necessary to do any real operations for updating 
}
