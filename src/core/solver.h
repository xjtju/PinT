#ifndef PinT_SOLVER_H_
#define PinT_SOLVER_H_ 

#include "PinT.h"
#include "Grid.h" 

/**
 * the abstract interface of all solvers for PinT framework
 */
class Solver {
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
    
    int inner_size; // guard cells not included 
    int outer_size; // guard cells included
    int &size = outer_size;      // alias of outer_size 

    int steps;     // the number of time steps in one time slice

public:
    Solver(PinT *conf, Grid *g){
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
    }

    // the template algorithm for linear system etc. 
    virtual void solve() = 0;

    // integrate over one time slice , the default implementation
    void evolve();
};

#endif
