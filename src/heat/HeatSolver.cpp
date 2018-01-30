#include "HeatSolver.h"

HeatSolver::HeatSolver(PinT *c, Grid *g) : DefaultSolver(c,g){ 
    setup();
}

HeatSolver::HeatSolver(PinT *c, Grid *g, bool isFS) : DefaultSolver(c,g,isFS){
    setup();
}

// set diffuse coefficient and tune the default parameter, problem specific
void HeatSolver::setup(){
    if(ndim==1) k = 0.061644; 
    if(ndim==2) k = 0.061644;
    if(ndim==3) k = 0.061644;
}

