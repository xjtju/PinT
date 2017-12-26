// PBiCBSTAB 
#include "PBiCGStab.h"

void PBiCGStab::solve(double* x, const int itmax){


    double rho0, rho, alpha, omega, beta;
    double res = 10.0;
   
    rho0 = alpha = omega = 1.0; 	
      
    clear_mem(v, size);
    clear_mem(p, size);

    cg_b(x);
    cg_rk(r, x, b); //init residual r = b - Ax 
	// choose an aribtrary vector r0_ such that (r0_, r) /=0, e.g. r0_ = r
	blas_cp(r0_, r, size);

    double b_nrm2 = blas_vdot(b, b, size);

    double tmp, tmp1, tmp2;	
     
    for(int i=1; i<=itmax; i++){
		rho = blas_vdot(r0_, r, size);
        beta = (rho/rho0) * (alpha/omega);
        cg_direct(p, r, v, beta, omega); 

	    preconditioner(p_,p); //solve Mp_=p, in some algorithm description, p_ is also denoted by y
	    cg_Xv(v,p_);        // v=Ap_	
        
        tmp = blas_vdot(r0_, v, size);
        alpha = rho / tmp; 
        // h = x + alpha*p_, is ||h|| is sufficiently small, can do...  
        
        blas_avpy(s, -alpha, v, r, size); //s = -av + r  or s = r - av
        // solve Ms_=s , in some algorithm description, s_ is also denoted by z
	    preconditioner(s_,s); 
        cg_Xv(t,s_); // t=Az
        // in preconditioner t = 1/M_1*t and s = 1/M_1*s
	    tmp1 = blas_vdot(t, s, size); // t dot z	
        tmp2 = blas_vdot(t, t, size); // t dot t 
		omega = tmp1 / tmp2;  // 
        
	    cg_xi(x, alpha, p_, omega, s_); // xi = xi_1 + alpha*y + omega*z; 	
	    blas_avpy(r, -omega, t, s, size); //r = s - omega*t	
         
		// ||r|| 
		tmp = blas_vdot(r, r, size);
        res = sqrt(tmp/b_nrm2); 
		if(res < eps) {
            //printf("%4d\n", i);
            break;
        }
        rho0 = rho;
	}
    
}

//p = r + beta * ( p - omg * v )
void PBiCGStab::cg_direct(double* p, double* r, double* v, double beta, double omega){
    for(int i=0; i<size; i++) {
        p[i] = r[i] + beta*(p[i] - omega*v[i]);
    }
}

// x = x + ay + bz  
void PBiCGStab::cg_xi(double* x, double alpha, double* y, double omega, double* z){
    for(int i=0; i<size; i++) {
        x[i] = x[i] + alpha*y[i] + omega*z[i];
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
