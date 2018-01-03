// PBiCBSTAB 
#include "PBiCGStab.h"

void PBiCGStab::solve(){

    double *x = fetch();  // get the grid variables

    double rho0, rho, alpha, omega, beta;
    double res = 1000.0;
   
    rho0 = alpha = omega = 1.0;     
      
    clear_mem(v, size);
    clear_mem(p, size);
    clear_mem(b, size);
    clear_mem(r, size);

    grid->guardcell(x);     
    cg_b(x);
    //grid->guardcell(b);
    cg_rk(r, x, b); //init residual r = b - Ax 
    //grid->guardcell(r);

    // choose an aribtrary vector r0_ such that (r0_, r) /=0, e.g. r0_ = r
    blas_cp(r0_, r, size);

    double tmp, tmp1, tmp2; 
    tmp = tmp1 = tmp2 = 0.0;  

    double b_nrm2 = blas_vdot(b, b, nx, nguard);
    grid->sp_allreduce(&b_nrm2);

    for(int i=1; i<=itmax; i++){
        rho = blas_vdot(r0_, r, nx, nguard);
        //printf("%d rho : %f\n",i,  rho);
        grid->sp_allreduce(&rho);

        beta = (rho/rho0) * (alpha/omega);
        cg_direct(p, r, v, beta, omega); 

        grid->guardcell(p);
        preconditioner(p_,p); //solve Mp_=p, in some algorithm description, p_ is also denoted by y

        grid->guardcell(p_);
        cg_Xv(v,p_);        // v=Ap_    
        
        tmp = blas_vdot(r0_, v, nx, nguard);
        grid->sp_allreduce(&tmp);
        alpha = rho / tmp; 
        // h = x + alpha*p_, is ||h|| is sufficiently small, can do...  
        
        blas_avpy(s, -alpha, v, r, nx, nguard); //s = -av + r  or s = r - av
        grid->guardcell(s);
        // solve Ms_=s , in some algorithm description, s_ is also denoted by z
        preconditioner(s_,s); 

        grid->guardcell(s_);
        cg_Xv(t,s_); // t=Az
        // in preconditioner t = 1/M_1*t and s = 1/M_1*s
        tmp1 = blas_vdot(t, s, size); // t dot z    
        grid->sp_allreduce(&tmp1);
        tmp2 = blas_vdot(t, t, size); // t dot t 
        grid->sp_allreduce(&tmp2);
        omega = tmp1 / tmp2;   
        
        cg_xi(x, alpha, p_, omega, s_); // xi = xi_1 + alpha*y + omega*z;   
        blas_avpy(r, -omega, t, s, size); //r = s - omega*t 
         
        // ||r|| 
        tmp = blas_vdot(r, r, nx, nguard);
        grid->sp_allreduce(&tmp);   
        res = sqrt(tmp/b_nrm2);

        //printf("PBI %d, %f\n",i, res);
        if(res < eps) {
            break;
        }
        rho0 = rho;
        grid->bc(x);
    }

    grid->guardcell(x);
    update(); // update grid variables
}

//
// vector operations can be accelerated by SIMD. 
//

//p = r + beta * ( p - omg * v )
void PBiCGStab::cg_direct1d(double* p, double* r, double* v, double beta, double omega){
    long idx;
    for(int i=nguard; i<nx+nguard; i++) {
        idx = grid->getInnerIdx(i);
        p[i] = r[idx] + beta*(p[i] - omega*v[idx]);
    }
}

// p : outer_size; r, v : inner_size
void PBiCGStab::cg_direct2d(double* p, double* r, double* v, double beta, double omega){
    long idx;
    for(int j=nguard; j<ny+nguard; j++)
        for(int i=nguard; i<nx+nguard; i++) {
            idx = grid->getInnerIdx(i,j); 
            p[i] = r[idx] + beta*(p[i] - omega*v[idx]);
    }
}

// x = x + ay + bz (x p_ s, outer_size) 
void PBiCGStab::cg_xi1d(double* x, double alpha, double* y, double omega, double* z){
    for(int i=nguard; i<nx+nguard; i++) {
        x[i] = x[i] + alpha*y[i] + omega*z[i];
    }
}

void PBiCGStab::cg_xi2d(double* x, double alpha, double* y, double omega, double* z){
    long idx;
    for(int j=nguard; j<ny+nguard; j++) 
        for(int i=nguard; i<nx+nguard; i++) {
            idx = j*sy + i;
            x[idx] = x[idx] + alpha*y[idx] + omega*z[idx];
    }
}

/**
  preconditioner, in current version , there is no precondition
  solve Mp_=p
  M = M1*M2 ~ A, there is a simple method for decomposing the matrix A when A has some regular pattern, 
  for example in HEAT diffusion:
  if we set  M1=I, M2=A, where I is Identity matrix. 
  than the preconditioner can be set to a light-weight but less precise solver than BiCG, and solving the same linear system: Ax=b, 
  moreover the following calculations of BiCG will become : 
    ==  reverse(M1)*t=I*t=t and reverse(M1)*s=I*s=s ==  
  so the matrix reverse and matrix multiplying vector calculations can also be skipped.      
*/
void PBiCGStab::preconditioner(double* p_, double* p){
    blas_cp(p_, p, size);
}

void PBiCGStab::update() {
    //in current implementation, the solver directly modify the grid variables, 
    //so it is not necessary to do any real operations for updating 
}
