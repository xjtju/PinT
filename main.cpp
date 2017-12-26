#include <mpi.h>

#include "common.h"
#include "HeatCG.h"
#include "HeatModel.h"

// integrate the target equation along one time slice,   
int integrator(PBiCGStab *solver,Model *model);

int main(int argc, char* argv[]) {
    
    int myid, numprocs;

    int spnum;    // space parallel 
    int tsnum;    // time slice   

    int c_steps;  // number of time steps of coarse solver within one time slice  
    int f_steps;  // fine solver 
    int c_dt;       // time step with of coarse solver
    int f_dt;       // fine solver

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    spnum = 1;
    tsnum = numprocs / spnum; 
    c_steps = 10;
    f_steps = 100;

    
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
