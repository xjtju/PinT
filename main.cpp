#include <mpi.h>

#include "common.h"
#include "PinT.h"
#include "HeatGrid.h"
#include "HeatSolver.h"

// the PinT algorithm template
int evolve(PinT *conf, Grid *fg, PBiCGStab *fs, Grid *gg, PBiCGStab *gf);

// integrate the target equation along one time slice,   
void integrate(Grid *g, PBiCGStab *solver, int steps);

int main(int argc, char* argv[]) {
    
    int myid, numprocs;

    int spnum;    // space parallel 
    int tsnum;    // time slice   
    
    MPI_Comm sp_comm; //space within the same time slice

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    PinT *conf = new PinT();
    conf->init();
    spnum = conf->subgrids;
    tsnum = conf->slices; 
    
    //check the configuration is consist with the real run time
    assert(tsnum == (numprocs/spnum));
    
    if(myid == 0) {
       conf->print(); 
    }
   
    int key=myid%spnum, color=myid/spnum;  
    int spid=0;
    
    //communicator for spatil parallel
    MPI_Comm_split(MPI_COMM_WORLD,color,key,&sp_comm);
    MPI_Comm_rank(sp_comm, &spid);

    // create Coarse/fine Grid and Solver
    HeatGrid *fg = new HeatGrid(
            conf->sub_nx, // grid size
            1,   //guard cells
            conf->dx, //cell size
            conf->f_dt //time step
            );
    fg->init();
    HeatSolver *fslv = new HeatSolver(fg);    
    fslv->set_eps(1.0e-6);
    fslv->set_itmax(10);

    HeatGrid *gg = new HeatGrid(
            conf->sub_nx,
            1,
            conf->dx,
            conf->c_dt
            );
    gg->init();
    HeatSolver *gslv = new HeatSolver(gg);
    gslv->set_eps(1.0e-6);
    gslv->set_itmax(10);

    int source, dest, tag;
    int ierr;
    int size = gg->size; 
    MPI_Request req;
    MPI_Status  stat;

    double *c_prev = alloc_mem(gg->size); // the structure holder for the coarse solver at the previous iteration of the the same slice 
    double *u_start= alloc_mem(gg->size); // the latest solution of the current time slice start point or the previous slice end 
    double *u_end  = alloc_mem(gg->size); // the latest solution of the current time slice end point or the next slice start 

    blas_cp(u_start, gg->x, size);
    if(myid/spnum >= 1) {
        source = myid - spnum;
        tag = myid*100;
        MPI_Recv(u_start, size, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &stat);
        // besides the first time slice, all others need to receive U^{0}_{n-1} as its start value  
    }
    
    blas_cp(gg->x, u_start, size); 
    //coarse 
    integrate(gg, gslv, conf->c_steps);
      
    if(myid/spnum<(tsnum-1)){
        dest = myid + spnum;
        tag  = (myid + spnum)*100;
        //send to next time slice
        ierr = MPI_Rsend(gg->x, size, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);    
    }
    blas_cp(u_end, gg->x, size); 

    int kpar=1; 
    double res_loc, res_sp, max_res;
    res_loc = res_sp = max_res = 0.0;
    for(; kpar<=conf->kpar_limit; kpar++)
    {
        // step1: fine solver parallel run based on U^{k-1}_{n-1} 
        blas_cp(fg->x, u_start, size);
        integrate(fg, fslv, conf->f_steps);
        // step2:
        blas_cp(c_prev, gg->x, size); 
        // step3:
        if(myid/spnum>=1){
	        source = myid - spnum;
	        tag    = myid*100 + kpar;
	        MPI_Recv(u_start, size, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &stat);
        } 
        // step4:
        blas_cp(gg->x, u_start, size); 
        integrate(gg, gslv, conf->c_steps);
        // step5: 
        blas_pint_sum(u_end, fg->x, gg->x, c_prev, &res_loc, size);  
        // step6: 
        if(myid/spnum< (tsnum-1)){
	        dest = myid + spnum;
	        tag  = (myid + spnum)*100 + kpar;
	        ierr = MPI_Send(u_end, size, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);    
        }
        
        //STEP7 gather residual
        MPI_Allreduce(&res_loc, &res_sp,  1, MPI_DOUBLE, MPI_SUM, sp_comm);
        MPI_Allreduce(&res_sp,  &max_res, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        
        if(myid==numprocs-1){
	        printf("id:%d, kpar:%d, res:%13.6e \n", myid,kpar,sqrt(max_res));
        }
     
        //STEP8
        if( sqrt(max_res) < conf->converge_eps) 
            break;   
    }
    
    if(myid==numprocs-1) { 
        for(int i=0; i<size; i++){
            printf("%f\n", u_end[i]);
        }
        printf("kpar=%d\n",kpar);
    }

    delete fslv;
    delete fg;
    delete gg;
    delete gslv;

    MPI_Barrier(MPI_COMM_WORLD); 

    return 0;
}

void integrate(Grid *g, PBiCGStab *solver, int steps){

    for(int i=0; i<steps; i++){
        g->bc();
        solver->solve();
    }
}

int evolve(PinT *conf, Grid *fg, PBiCGStab *fs, Grid *gg, PBiCGStab *gf)
{
    for(int i=0; i<conf->Nt; i++){
        fg->bc();
        fs->solve();
    }

    return 0;
}
