#include "common.h"
#include "HeatCG.h"
#include "HeatModel.h"

// integrate the target equation along one time slice,   
int integrator(PBiCGStab *solver,Model *model);

int main(int argc, char* argv[]) {

    int size = NX + 2*NGuard;
    Model *model = new HeatModel(NX, NGuard);
    PBiCGStab *solver = new HeatCG(NX, NGuard, EPS);    

    model->init_x();
    
    integrator(solver, model); 
    
    
    
    double *x = model->x;
    for(int i=NGuard; i<=NX; i++){
        printf("%f\n", x[i]);
    }

    delete solver;
    delete model;

    return 0;
}


int integrator(PBiCGStab *solver, Model *model){

    for(int i=0; i<NT; i++){
        model->bc();
        solver->solve(model->x,10);
    }
    return 0;
}
