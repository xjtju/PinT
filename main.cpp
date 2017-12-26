#include <mpi.h>

#include "common.h"
#include "PinT.h"
#include "HeatGrid.h"
#include "HeatCG.h"

// integrate the target equation along one time slice,   
int integrator(PBiCGStab *solver,PinT *conf, Grid *grid);

int main(int argc, char* argv[]) {
    
    int myid, numprocs;

    int spnum;    // space parallel 
    int tsnum;    // time slice   

    int c_steps;  // number of time steps of coarse solver within one time slice  
    int f_steps;  // fine solver 
    int c_dt;       // time step with of coarse solver
    int f_dt;       // fine solver
    /*
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    spnum = 1;
    tsnum = numprocs / spnum; 
    */
    PinT *conf = new PinT();
    conf->init();
    
    HeatGrid *fg = new HeatGrid(
            conf->sub_nx,
            1,
            conf->dx,
            conf->f_dt
            );

    fg->init();
    int nx = fg->nx;
    int size = fg->size; 
    int nguard = fg->nguard;

    PBiCGStab *solver = new HeatCG(fg, EPS);    


    integrator(solver, conf, fg); 
    
    
    
    double *x = fg->x;
    for(int i=nguard; i<=nx; i++){
        printf("%f\n", x[i]);
    }

    delete solver;
    delete fg;

    return 0;
}


int integrator(PBiCGStab *solver, PinT *conf, Grid *grid){

    for(int i=0; i<conf->Nt; i++){
        grid->bc();
        solver->solve(grid->x,10);
    }
    return 0;
}
