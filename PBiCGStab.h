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
#include "Grid.h"

class PBiCGStab {

protected:
    int nx; 
    int nguard;
    int size;

    double eps = 1.0e-6;
    int itmax = 10; 

    double *r0_, *r;  // 0 ... itmax  
    double *p, *p_;   // conjugated direction
	double *v;        //direction
    double *s, *s_, *t;
    double *b;       //RHS

    Grid *grid;
    
public:

    inline void set_eps(double eps) { this->eps = eps; }
    inline void set_itmax(int iter) { this->itmax = iter;}

    PBiCGStab(Grid *g){
        this->grid = g;
        this->nx = g->nx;
        this->nguard = g->nguard;
        this->size = g->size; 
        this->eps = eps;

        r0_ = alloc_mem(size);
        r   = alloc_mem(size);
        p   = alloc_mem(size);
        p_  = alloc_mem(size);
        v   = alloc_mem(size);
        s   = alloc_mem(size);
        s_  = alloc_mem(size);
        t   = alloc_mem(size);

        b   = alloc_mem(size);
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


    // the template algorithm of PBiCBSTAB
    void solve();

    // x = x + ay + bz  
    void cg_xi(double *x, double alpha, double *y, double omega, double *z);

    //p = r + beta * ( p - omg * v )
    void cg_direct(double* p, double* r, double* v, double beta, double omega); 

    // preconditioner: solve Mp_=p
    void preconditioner(double *p_, double *p);
    
    double* fetch(); // fetch physical variables into solver
    void update();  // update physical variables in grid 

    /***** stencil related functions *****/
    //calcaluate the residual r = b - Ax 
    virtual void cg_rk(double *r, double *x, double *b) = 0;
    // matrix * vector, v = A*y 
    virtual void cg_Xv(double *v, double *y) = 0; 	
    // set the b in Ax=b 
    virtual void cg_b(double *x)=0;
};

#endif
