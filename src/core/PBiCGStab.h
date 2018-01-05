#ifndef _PinT_PBiCBStab_H_
#define _PinT_PBiCBStab_H_ 1

/**
 * the template of Biconjugate gradient stabilized method for solving Ax=B.
 * Except the unknown vector x, all the temporary variables are managed by the template including the RHS b.
 *
 * Stencil related routines are declared by virtual functions, the concrete class must implement them  
 * 
 * For easy readning, most functions pass parameter explicitly though it is not necessary because all variables are shared within class
 */

#include "common.h"
#include "blas.h"
#include "blascg.h"
#include "Solver.h"

class PBiCGStab : public Solver{

protected:

    double eps = 1.0e-6;
    int itmax = 20; 

    double *r0_, *r ;  // 0 ... itmax  
    double *p, *p_ ;   // conjugated direction, out size
    double *v;         //direction
    double *s, *s_ ;   // out size 
    double *t;
    double *b;         //RHS

public:

    inline void set_eps(double eps) { this->eps = eps; }
    inline void set_itmax(int iter) { this->itmax = iter;}

    PBiCGStab(PinT* c, Grid *g):Solver(c,g){
        b   = alloc_mem(size);

        r0_ = alloc_mem(size);
        r   = alloc_mem(size);
        v   = alloc_mem(size);
        t   = alloc_mem(size);

        p   = alloc_mem(size);
        p_  = alloc_mem(size);
        s   = alloc_mem(size);
        s_  = alloc_mem(size);
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

       free_mem(b);
    }


    // preconditioner: solve Mp_=p
    void preconditioner(double *p_, double *p);
    
    virtual double* fetch() = 0; // fetch physical variables into solver
    void update() ;  // update physical variables in grid, reserved blank functions. 
    
    // the template algorithm of PBiCBSTAB
    void solve();


    /***** stencil related functions *****/
    //calcaluate the residual r = b - Ax 
    inline void cg_rk(double *r, double *x, double *b){
        if(ndim==1) cg_rk1d(r, x, b);
        else if(ndim==2) cg_rk2d(r, x, b);
    }
    virtual void cg_rk1d(double *r, double *x, double *b)=0;
    virtual void cg_rk2d(double *r, double *x, double *b)=0;

    // matrix * vector, v = A*y 
    inline void cg_Xv(double *v, double *y) {
        if(ndim==1) cg_Xv1d(v, y);
        else if(ndim==2) cg_Xv2d(v, y);
    }
    virtual void cg_Xv1d(double *v, double *y)=0;
    virtual void cg_Xv2d(double *v, double *y)=0;

    // set the b in Ax=b, RHS 
    inline void cg_b(double *x) {
        if(ndim==1) cg_b1d(x);
        else if(ndim==2) cg_b2d(x);
    }
    virtual void cg_b1d(double *x)=0;
    virtual void cg_b2d(double *x)=0;

   
    /**** BLAS related functions *****/
    // vector production 
    inline double cg_vdot(double *t, double *z){
        double tmp;
        if(ndim==1) tmp = blas_vdot(t, z, nx, nguard); 
        else if(ndim==2) 
            blas_vdot_2d_(grid->nxyz, &nguard, t, z, &tmp); 
        return tmp;
    }

    // s = r - alpha*v or r = s - omega*t 
    inline void cg_avpy(double* s, double alpha, double* r, double* v) {
        if(ndim==1) blas_avpy(s, -alpha, v, r, nx, nguard);
        else if (ndim==2) blas_avpy_2d_(grid->nxyz, &nguard, &alpha, s, r, v); 
    }

    /**** CG special functions *****/
    // in practice, from 2D, fortran is used.
    //
    // x = x + ay + bz  
    inline void cg_xi(double *x, double alpha, double *y, double omega, double *z){
        if(ndim==1) 
            cg_xi1d(x, y, z, alpha, omega);
        else if(ndim==2) 
            cg_xi2d_(grid->nxyz, &nguard, x, y, z, &alpha, &omega);
    }
    void cg_xi1d(double *x, double *y, double *z, double alpha, double omega);
    void cg_xi2d(double *x, double *y, double *z, double alpha, double omega);

    //p = r + beta * ( p - omg * v )
    inline void cg_direct(double* p, double* r, double* v, double beta, double omega){
        if(ndim==1)
            cg_direct1d(p, r, v, beta, omega);
        else if(ndim==2)
            cg_direct2d_(grid->nxyz, &nguard, p, r, v, &beta, &omega);
    }
    void cg_direct1d(double* p, double* r, double* v, double beta, double omega);
    void cg_direct2d(double* p, double* r, double* v, double beta, double omega);

};

#endif
