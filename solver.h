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

    int nx; 
    int nguard;
    int size;
    int steps;

public:
    Solver(PinT *conf, Grid *g){
        this->conf = conf; 
        this->grid = g;

        this->nx = g->nx;
        this->nguard = g->nguard;
        this->size = g->size;
    }

    // the template algorithm for linear system etc. 
    virtual void solve() = 0;

    // integrate over one time slice , the default implementation
    void evolve(){
        for(int i=0; i<steps; i++){
            grid->bc();
            solve();
        }
    }
};

#endif
