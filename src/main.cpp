#include <mpi.h>

#include "common.h"
#include "Driver.h"
#include "Solver.h"
#include "HeatGrid.h"
#include "HeatSolverF.h"

// the PinT algorithm template
int evolve(PinT *conf, Grid *fg, PBiCGStab *fs, Grid *gg, PBiCGStab *gf);

/**
 * the ini configuration file as the first command line parameter
 */
int main(int argc, char* argv[]) {
    
    char* ini_file = (char*)"pint.ini";  // default ini file
   
    if(argc>1) {
        ini_file = argv[1];      
    }
   
    int myid, numprocs;

    int spnum;    // space parallel 
    int tsnum;    // time slice   
    
    MPI_Comm sp_comm; //space within the same time slice

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    Driver driver;

    //load global configuration 
    PinT *conf = PinT::instance();
    driver.init(ini_file,conf);

    spnum = conf->spnum;
    tsnum = conf->tsnum; 
    
    //check the configuration is consist with the real run time
    if(tsnum != (numprocs/spnum)) {
        driver.Abort("configuration inconsist : %s. \n" "(the space num)*(the time num) should equal with the total process num");
    }

    int key=myid%spnum, color=myid/spnum;  
    int spid=0;
    //communicator for spatil parallel
    MPI_Comm_split(MPI_COMM_WORLD,color,key,&sp_comm);
    MPI_Comm_rank(sp_comm, &spid);

    // create the grid/mesh and solver 
    Grid *g = new HeatGrid(conf);
    g->init();

    Solver *F = new HeatSolverF(conf,g);    
    Solver *G = new HeatSolverC(conf,g);
    // for convience only, set the pointer to grid inner variables
    double *u_cprev = g->u_cprev;  
    double *u_start= g->u_start; 
    double *u_end  = g->u_end;  
    double *u_c    = g->u_c;
    double *u_f    = g->u_f;

    int source, dest, tag;
    int ierr;
    int size = g->size; 
    MPI_Request req;
    MPI_Status  stat;

    blas_cp(u_start, u_c, size);
    if(myid/spnum >= 1) {
        source = myid - spnum;
        tag = myid*100;
        MPI_Recv(u_start, size, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &stat);
        // besides the first time slice, all others need to receive U^{0}_{n-1} as its start value  
    }
    
    g->bc();
    blas_cp(u_c, u_start, size); 
    //coarse 
    G->evolve();
      
    g->bc();
    if(myid/spnum<(tsnum-1)){
        dest = myid + spnum;
        tag  = (myid + spnum)*100;
        //send to next time slice
        ierr = MPI_Rsend(u_c, size, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);    
    }
    blas_cp(u_end, u_c, size); 

    int kpar=1; 
    double res_loc, res_sp, max_res;
    res_loc = res_sp = max_res = 0.0;

    for(; kpar<=conf->kpar_limit; kpar++)
    {
        res_loc = 0.0;
        res_sp  = 0.0;
        max_res = 0.0;
        g->bc();
        // step1: fine solver parallel run based on U^{k-1}_{n-1} 
        blas_cp(u_f, u_start, size);
        F->evolve();
        g->bc();
        // step2:
        blas_cp(u_cprev, u_c, size); 
        // step3:
        if(myid/spnum>=1){
	        source = myid - spnum;
	        tag    = myid*100 + kpar;
	        MPI_Recv(u_start, size, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &stat);
        } 
        // step4:
        blas_cp(u_c, u_start, size); 
        G->evolve();
        g->bc();
        // step5: 
        blas_pint_sum(u_end, u_f, u_c, u_cprev, &res_loc, size);  
        // step6: 
        if(myid/spnum< (tsnum-1)){
	        dest = myid + spnum;
	        tag  = (myid + spnum)*100 + kpar;
	        ierr = MPI_Send(u_end, size, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);    
        }
        
        //STEP7 gather residual
        MPI_Allreduce(&res_loc, &res_sp,  1, MPI_DOUBLE, MPI_SUM, sp_comm);
        MPI_Allreduce(&res_sp,  &max_res, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        
        //res_loc = sqrt(res_loc);
        //if(myid==0){
	    //    printf("kpar:%d, myid:%d, loc_res:%13.6e \n\n", kpar, myid, res_loc);
        //}
	    //printf("kpar:%d, myid:%d, res_loc:%13.6e \n", kpar, myid, res_loc);

        max_res = sqrt(max_res/tsnum);
        if(myid==numprocs-1){
	        printf("kpar:%d, myid:%d, max_res:%13.6e \n\n", kpar, myid, max_res);
        }
     
        //STEP8
        if( max_res < conf->converge_eps) 
            break;   
    }
    
    if(myid==numprocs-1) { 
        for(int i=0; i<size; i++){
            printf("%f\n", u_end[i]);
        }
        printf("kpar=%d\n",kpar);
    }

    delete F;
    delete G;
    delete g;

    MPI_Barrier(MPI_COMM_WORLD);

    return 0;
}

