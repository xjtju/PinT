#include "Driver.h"
#include "blas.h"
/**
 * load parammeter configuration file and check their consistency 
 */
void Driver::init(int argc, char* argv[]){

    char* ini_file = (char*)"pint.ini";  // default ini file
   
    if(argc>1) {
        ini_file = argv[1];  // the ini configuration file as the first command line parameter
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   
    conf = PinT::instance();

    INFO("Initialization file is  %s \n", ini_file);

    if (ini_parse(ini_file, handler, conf) < 0) {
        Abort("Can't load ini file : %s .\n", ini_file);
    } else {
        conf->init();
        if(myid==0) conf->print();
    }

    spnum = conf->spnum;
    tsnum = conf->tsnum;
    mytid = myid / spnum;

    //check the configuration is consist with the real run time
    if(tsnum != (numprocs/spnum)) {
        Abort("configuration inconsistent : %s. \n", "(the space num)*(the time num) should be equal with the total process num");
    }

    if(spnum != (conf->spnumx*conf->spnumy*conf->spnumz)) {
        Abort("configuration inconsistent : %s. \n", "(the space num) should be equal with value multiplied from all directions");
    }

    int key=myid%spnum, color=myid/spnum; 
    //communicator for spatil parallel
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &sp_comm);
    MPI_Comm_rank(sp_comm, &mysid);
    
    conf->sp_comm = &sp_comm;
    conf->myid = myid;
    conf->mytid = mytid;  
    conf->mysid = mysid;
}

// the PinT algorithm template
void Driver::evolve(Grid* g, Solver* G, Solver* F){
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

    blas_cp(u_start, u_c, size); // is not necessary, because all vectors have the same values at the init
     // except the first time slice, all others need to receive U^{0}_{n-1} as its start value  
    if(!isFirstSlice(myid)) {
        source = myid - spnum;
        tag = myid*100;
        MPI_Recv(u_start, size, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &stat);
        g->bc();
    }
    
    //coarse 
    blas_cp(u_c, u_start, size); 
    G->evolve();
    g->bc();

    // except the last time slice, all others need to send the coarse(estimate) value to its next slice  
    if(!isLastSlice(myid)){
        dest = myid + spnum;
        tag  = (myid + spnum)*100;
        //send to next time slice
        ierr = MPI_Rsend(u_c, size, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);    
    }
    
    //blas_cp(u_end, u_c, size); 

    double res_loc, res_sp, max_res;
    res_loc = res_sp = max_res = 0.0;

    for(int k=1; k<=conf->kpar_limit; k++)
    {
        res_loc = 0.0;
        res_sp  = 0.0;
        max_res = 0.0;

        kpar = kpar + 1;
        g->bc();
        
        // step1:
        blas_cp(u_cprev, u_c, size); //this step is not necessary at the following of fine solver 

        // step2: fine solver parallel run based on U^{k-1}_{n-1}
        blas_cp(u_f, u_start, size);
        F->evolve();
        g->bc();

        if(kpar == 1) {
            blas_cp(u_end, u_f, size); 
        }
        g->bc();

        // step3:
        if(!isFirstSlice(myid)){
	        source = myid - spnum;
	        tag    = myid*100 + kpar;
	        MPI_Recv(u_start, size, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &stat);
            g->bc();
        } 
        // step4:
        blas_cp(u_c, u_start, size); 
        G->evolve();
        g->bc();
        
        // step5: 
        blas_pint_sum(u_end, u_f, u_c, u_cprev, &res_loc, size, true);  
        // step6: 
        if(!isLastSlice(myid)){
	        dest = myid + spnum;
	        tag  = (myid + spnum)*100 + kpar;
	        ierr = MPI_Send(u_end, size, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);    
        }
        
        //STEP7 gather residual
        MPI_Allreduce(&res_loc, &res_sp,  1, MPI_DOUBLE, MPI_SUM, sp_comm);
        MPI_Allreduce(&res_sp,  &max_res, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        
        monitorResidual(u_c,u_cprev,u_end,u_f,res_loc,max_res,size);

        max_res = sqrt(max_res/tsnum);
        //STEP8
        if( max_res < conf->converge_eps) 
            break;   
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

void Driver::finalize() {
    MPI_Finalize();
}

/**
 * Due to the "Machine Epsilon" or rounding error, residual calculation must be very careful.
 * Sometimes, in theory, the residual should be ZERO, but in practice, the calculation value is not ZERO, 
 * despite it is very small, it will has an unignorable impact on convergency due to  accumulating effect.  
 *
 * the fine solver starts from the first time slice, as it evolving, for all the slices which it has walked, 
 * it is no need to calcaluate again after all next processes.  
 * That is the accumulative residual before the time slice should be ZERO according PARAREAL's theory..  
 */
void Driver::monitorResidual(double* u_c, double* u_cprev, double* u_end, double* u_f, double res_loc, double max_res,int size ){

    bool debug = false;
    if(debug && (myid==0)){
        res_loc = sqrt(res_loc);
        double cdist = blas_vdist(u_c,u_cprev,size); 
        double fdist = blas_vdist(u_end,u_f,size); 
        printf("kpar:%d, myid:%d, cvdist:%13.8e, fvdist:%13.8e , loc_res:%13.8e\n", kpar, myid, cdist, fdist, res_loc); 
    }

    if(mytid == kpar-1){
        if(res_loc > 0 ) {
            res_loc = sqrt(res_loc);
            double cdist = blas_vdist(u_c,u_cprev,size); 
            double fdist = blas_vdist(u_end,u_f,size); 
            WARN("Local residual should be ZERO, but it's NOT, details : \n \t kpar:%d, myid:%d, cvdist:%13.8e, fvdist:%13.8e , loc_res:%13.8e . \n\n", kpar, myid, cdist, fdist, res_loc); 
        }
    }

    if(myid==numprocs-1){
	    printf("kpar:%d, myid:%d, max_res:%13.8e \n", kpar, myid, max_res);
    }
}

void Driver::INFO(const char* fmt, ...) {

    if(myid != 0) return ;
     
    int ret;
    va_list args;

    va_start(args, fmt);
    fputs("INFO : ", stdout);
    ret = vfprintf(stdout, fmt, args);
    va_end(args);
}

void Driver::WARN(const char* fmt, ...) {

    if(myid != 0) return ;
     
    int ret;
    va_list args;

    va_start(args, fmt);
    fputs("WARN : ", stderr);
    ret = vfprintf(stderr, fmt, args);
    va_end(args);
}


void Driver::Abort(const char* fmt, ...) {

    va_list args;
    
    va_start(args, fmt);
    fputs("ERROR: ", stderr);
    vfprintf(stderr, fmt, args);
    va_end(args);

    MPI_Finalize();

    exit(1);
}

