#include "Solver.h"

void Solver::evolve(){
    for(int i=0; i<steps; i++){
//            grid->bc();
        solve();
    }
}
