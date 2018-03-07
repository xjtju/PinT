#ifndef PinT_LS_H_
#define PinT_LS_H_  1

#include "PinT.h"
#include "Grid.h" 

/**
 * the abstract linear solver for PinT framework
 * 
 * two algorithms are provided in current version
 *   0. BiCG
 *   1. SOR
 *
 *  NOTE:
 *    I had considered creating an independent section for each solver,
 *    but it will introduce confusion if multiple instances of a solver need different settings.   
 */

class LS {

protected:
    Grid *grid;
    PinT *conf;
    
    // all the following variables are the same with the corresponding one in the Grid
    // holding the most being used variables here is just for convenience only 
    int ndim;

    int nx;
    int ny;
    int nz;

    int sx;
    int sy;
    int sz;

    int nguard;
    
    size_t inner_size; // guard cells not included 
    size_t outer_size; // guard cells included
    size_t &size = outer_size;      // alias of outer_size 
    
    // these control parammeters can be over-writen by sub classes
    // the max iterations of iterative algorithms 
    int itmax = 20; 
    // the convergence error tolerance of iterative algorithms 
    double eps = 1.0e-6;

    int force_abort = 0; 

    // NOTE : 
    // Neither of preconditioner and omega is a common attribute for linear solver 
    // but it will be very convenient for debugging if adding them to the .INI file   
    bool   isPrecond  = false; // if preconditioner is needed.
    double relaxfactor= 1.5;     // For SOR only 

public:
    LS(PinT *conf, Grid *g){
        this->conf = conf; 
        this->grid = g;

        this->ndim = g->ndim;

        this->nx = g->nx;
        this->ny = g->ny;
        this->nz = g->nz;

        this->sx = g->sx;
        this->sy = g->sy;
        this->sz = g->sz;

        this->nguard = g->nguard;

        this->inner_size = g->inner_size; 
        this->outer_size = g->size;

        this->itmax = conf->ls_itmax;
        this->eps = conf->ls_eps;
        this->force_abort = conf->ls_abort;
        this->isPrecond  = conf->ls_precond; 
        this->relaxfactor = conf->ls_relaxfactor;

    }

    void set_eps(double eps) { this->eps = eps; }
    void set_itmax(int iter) { this->itmax = iter;}
    void set_precond(bool flag) { this->isPrecond = flag;}

    virtual ~LS(){ }
    // the template algorithm for linear system etc. 
    // return the iteration count 
    virtual int solve(double *x, double *b, double *bcp) {
        if(conf->myid == 0)
           fprintf(stderr, "WARN : blank solve function is used, NOTHING will be updated!\n" );
    }  
};

#endif
