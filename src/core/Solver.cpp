#include "Solver.h"


/**
 * for one time slice iteration (Kpar)
 */
void Solver::evolve(){

    //steps = 1; // for DEBUG 
    for(int i=0; i<steps; i++){
        grid->bc();
        solve();
    }
}
