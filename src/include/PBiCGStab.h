#ifndef _PinT_PBiCBStab_H_
#define _PinT_PBiCBStab_H_ 1

/**
 * the template of Biconjugate gradient stabilized method for solving Ax=B.
 * Except the unknown vector x and the RHS b, all the temporary variables are managed by the template.
 *
 * Stencil matrix must be provided by the caller.
 * 
 * For easy reading, most functions pass parameter explicitly, though it is not necessary because all variables are shared within class
 */

#include "common.h"
#include "blas.h"
#include "blascg.h"
#include "LS.h"
#include "Driver.h"

class PBiCGStab : public LS{

private:
    void init(); // init the BiCG-specific variables

protected:
    double *r0_, *r ;  // 0 ... itmax  
    double *p, *p_ ;   // conjugated direction, out size
    double *v;         //direction
    double *s, *s_ ;   // out size 
    double *t;
    //double *b;         //RHS

    // these control parammeters can be over-writen by sub classes
    double eps = 1.0e-6;
    int itmax = 20; 
    double sor_omg = 1.7;     // SOR : relaxation factor, used for preconditioner, 
    bool   isPrecond = false; // in the current version, no preconditioner is GOOD.
    
public:

    inline void set_eps(double eps) { this->eps = eps; }
    inline void set_itmax(int iter) { this->itmax = iter;}

    PBiCGStab(PinT* c, Grid *g):LS(c,g){
        init();
    }
    
    ~PBiCGStab() {
       free_mem(r0_);
       free_mem(r);
       free_mem(p);
       free_mem(p_);
       free_mem(v);
       free_mem(s);
       free_mem(s_);
       free_mem(t);

       if(grid->myid==0 && conf->verbose)
       printf("INFO: the memory allocated by PBiCGStab solver has been released.\n");
    }


    // the template algorithm of PBiCBSTAB
    void solve(double *x, double *b, double *bcp);

    // preconditioner: solve Mp_=p, 
    // if sub class uses preconditioner, they must implement the virtual stencil function.  
    void preconditioner(double *p_, double *p, bool isPrecond);

    // the wrapper red-black successive over-relaxation (sor2_core)
    void sor2(double *p_, double *p, int lc_max, bool checkCnvg);

    inline void sor2_core(double *p_, double *p, int *color){
        /* 
        if(ndim==1) sor2_core_1d(p_, p, color);
        else if(ndim==2) sor2_core_2d(p_, p, color);
        else if(ndim==3) { printf("3D is not finished\n"); }; 
        */
    }
    /*
    virtual void sor2_core_1d(double *p_, double *p, int *color) = 0;
    virtual void sor2_core_2d(double *p_, double *p, int *color) = 0;
    */

    /***** stencil related functions, specific problem must provide the stencil matrix (bcp) *****/
    //calcaluate the residual r = b - Ax 
    inline void cg_rk(double *r, double *x, double *b, double *bcp){
        if(ndim==3)      cg_rk3d_(grid->nxyz, &nguard, r, x, b, bcp);
        else if(ndim==2) cg_rk2d_(grid->nxyz, &nguard, r, x, b, bcp);
        else if(ndim==1) cg_rk1d_(grid->nxyz, &nguard, r, x, b, bcp);
    }

    // matrix * vector, v = A*x 
    inline void cg_ax(double *v, double *x, double *bcp) {
        if(ndim==3)      cg_ax3d_(grid->nxyz, &nguard,v, x, bcp);
        else if(ndim==2) cg_ax2d_(grid->nxyz, &nguard,v, x, bcp);
        else if(ndim==1) cg_ax1d_(grid->nxyz, &nguard,v, x, bcp);
    }
   
    /**** BLAS related functions *****/
    // vector production 
    inline double cg_vdot(double *t, double *z){
        double val;
        if(ndim==1)      blas_dot_1d_(grid->nxyz, &nguard, t, z, &val ); 
        else if(ndim==2) blas_dot_2d_(grid->nxyz, &nguard, t, z, &val ); 
        else if(ndim==3) blas_dot_3d_(grid->nxyz, &nguard, t, z, &val ); 
        return val;
    }

    // s = r - alpha*v or r = s - omega*t 
    inline void cg_avpy(double* s, double alpha, double* r, double* v) {
        if(ndim==1)       blas_avpy_1d_(grid->nxyz, &nguard, &alpha, s, r, v); 
        else if (ndim==2) blas_avpy_2d_(grid->nxyz, &nguard, &alpha, s, r, v); 
        else if (ndim==3) blas_avpy_3d_(grid->nxyz, &nguard, &alpha, s, r, v); 
    }

    /**** CG special functions *****/
    // in practice, it is better to use fortran from 2D. 
    //
    // x = x + ay + bz  
    inline void cg_xi(double *x, double alpha, double *y, double omega, double *z){
        if(ndim==1)      cg_xi_1d_(grid->nxyz, &nguard, x, y, z, &alpha, &omega);
        else if(ndim==2) cg_xi_2d_(grid->nxyz, &nguard, x, y, z, &alpha, &omega);
        else if(ndim==3) cg_xi_3d_(grid->nxyz, &nguard, x, y, z, &alpha, &omega);
    }
    //p = r + beta * ( p - omg * v )
    inline void cg_direct(double* p, double* r, double* v, double beta, double omega){
        if(ndim==1)      cg_direct_1d_(grid->nxyz, &nguard, p, r, v, &beta, &omega);
        else if(ndim==2) cg_direct_2d_(grid->nxyz, &nguard, p, r, v, &beta, &omega);
        else if(ndim==3) cg_direct_3d_(grid->nxyz, &nguard, p, r, v, &beta, &omega);
    }

    // the following four methods are only used for debug
    void cg_xi1d(double *x, double *y, double *z, double alpha, double omega);
    void cg_xi2d(double *x, double *y, double *z, double alpha, double omega);
    void cg_direct1d(double* p, double* r, double* v, double beta, double omega);
    void cg_direct2d(double* p, double* r, double* v, double beta, double omega);
};

#endif
