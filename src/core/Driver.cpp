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
    if(argc>2){
        jobid = argv[2];  // the job id for profiling output 
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   
    conf = PinT::instance();

    INFO("Initialization file is  %s \n", ini_file);
    int ini_ret = conf->init(ini_file);
    if (ini_ret >= 0) {
        if(myid==0) conf->print();
    } else Abort("Can't load ini file or ini error : %s .\n", ini_file);

    conf->myid = myid;
    conf->numprocs = numprocs;
    spnum = conf->spnum;
    tsnum = conf->tsnum;
    smlr  = conf->smlr;
    mytid = myid / spnum;
    pipelined = conf->pipelined;
    
    //check the configuration is consist with the real run time
    bool checkOK = conf->check();
    if(!checkOK)
        Abort("configuration inconsistent : the program is forcely stopped due to the previous WARNs or ERRORs, please check the ini file and the run time parameters again!\n\n");

    int key=myid%spnum, color=myid/spnum; 
    //communicator for spatil parallel
    MPI_Comm_split(MPI_COMM_WORLD, color, key, &sp_comm);
    MPI_Comm_rank(sp_comm, &mysid);
    
    conf->sp_comm = &sp_comm;
    conf->myid = myid;
    conf->mytid = mytid;  
    conf->mysid = mysid;

    monitor.initialize();
    monitor.setRankInfo(myid);
}

// the PinT algorithm template
// NOTE:
//  the guardcell synchonization is the task of space parallel (Grid), not time parallel (Driver) 
//  according to the design principle of our framework. 
//  Guardcell synchonization is not necessary to be called explictly here except the end of the time loop    
//  for the consistency of output result with guardcell. 
void Driver::evolve(Grid* g, Solver* G, Solver* F){
    // output the initial values for debugging or ...
    if(conf->dump_init) {
        g->output_global_h5("ini");
        INFO("\ninitial data has been dumped to ...all.ini.h5 \n\n");
    }
    // for convience only, set the pointer to grid inner variables
    double *u_cprev = g->u_cprev;  
    double *u_start= g->u_start; 
    double *u_end  = g->u_end;  
    double *u_c    = g->u_c;  
    double *u_f    = g->u_f; 

    int source, dest, tag;
    int ierr;
    size_t size = g->size; 
    MPI_Request req;
    MPI_Status  stat;

    // except the first time slice, all others need to receive U^{0}_{n-1} as its start value  
    monitor.start(Monitor::RECV);
    if(!isFirstSlice()) {
        source = myid - spnum;
        tag = myid*100;
        MPI_Recv(u_start, size, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &stat);
    }
    monitor.stop(Monitor::RECV);

    //g->guardcell(u_start); // not necessary, the solver will do guardcell at the end of each time step   

    //coarse 
    monitor.start(Monitor::CSolver);
    blas_cp_(u_c, u_start, &size);  
    G->evolve();
    monitor.stop(Monitor::CSolver);
    
    monitor.start(Monitor::SEND);
    // except the last time slice, all others need to send the coarse(estimate) value to its next slice  
    if(!isLastSlice()){
        dest = myid + spnum;
        tag  = (myid + spnum)*100;
        //send to next time slice
        ierr = MPI_Rsend(u_c, size, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);    
    }
    monitor.stop(Monitor::SEND);

    // except the last time slice, all others need to send the coarse(estimate) value to its next slice  
    //blas_cp_(u_end, u_c, size); 

    double res_loc, res_sp, max_res;
    res_loc = res_sp = max_res = 0.0;
    for(int k=1; k<=conf->kpar_limit; k++)
    {
        res_loc = 0.0;
        res_sp  = 0.0;
        max_res = 0.0;

        kpar = kpar + 1;
        
        // step1:
        blas_cp_(u_cprev, u_c, &size); //this step is not necessary at the following of fine solver 

        // step2: fine solver parallel run based on U^{k-1}_{n-1}
        monitor.start(Monitor::FSolver);
        blas_cp_(u_f, u_start, &size);
        F->evolve();
        monitor.stop(Monitor::FSolver);

        if(kpar == 1) {
            blas_cp_(u_end, u_f, &size); 
        }

        // step3:
         monitor.start(Monitor::RECV); 
        if(!isFirstSlice()){
	        source = myid - spnum;
	        tag    = myid*100 + kpar;
	        MPI_Recv(u_start, size, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &stat);
        }
        monitor.stop(Monitor::RECV);
        // step4:
        
        monitor.start(Monitor::CSolver);
        blas_cp_(u_c, u_start, &size); 
        G->evolve();
        monitor.stop(Monitor::CSolver);
        
        // step5: 
        pint_sum(g, u_end, u_f, u_c, u_cprev, &res_loc, &smlr);  
        
        // step6:
        monitor.start(Monitor::SEND); 
        if(!isLastSlice()){
	        dest = myid + spnum;
	        tag  = (myid + spnum)*100 + kpar;
	        ierr = MPI_Send(u_end, size, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);    
        }
        monitor.stop(Monitor::SEND); 

        //STEP7 gather residual
        monitor.start(Monitor::RES);
        g->sp_allreduce(&res_loc, &res_sp);
        //MPI_Allreduce(&res_loc, &res_sp,  1, MPI_DOUBLE, MPI_SUM, sp_comm);
        if(pipelined == 0) {
            g->allreduce(&res_sp, &max_res, MPI_MAX);  // time space is enough ?!
            //MPI_Allreduce(&res_sp,  &max_res, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        }
        monitorResidual(g, res_loc, max_res, size);
        monitor.stop(Monitor::RES);

         //STEP8
        if(pipelined == 0){
            max_res = sqrt(max_res);
            if( max_res < conf->converge_eps) {
                if(myid==numprocs-1) printf("Parareal is converged at %d iteration, res = %e .\n", k , max_res);
                break;   
            }
        }
    }
    g->guardcell(u_end);
    MPI_Barrier(MPI_COMM_WORLD);
}

void Driver::finalize() {
    char fname[40];

    monitor.gather();

    memset(fname, 0, sizeof(char)*40);
    sprintf(fname, "%s_%s_basic.txt", conf->monitor_pre, jobid); 
    monitor.print(fname, "The TAO of Programming", "The PinT performance test framework");

    memset(fname, 0, sizeof(char)*40);
    sprintf(fname, "%s_%s_detail.txt", conf->monitor_pre, jobid); 
    monitor.printDetail(fname);

    monitor.print(stdout, "The TAO of Programming", "The PinT performance test framework");

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
void Driver::monitorResidual(Grid *g, double res_loc, double max_res,int size ){

    bool debug = false;
    double *u_c    = g->u_c;
    double *u_cprev = g->u_cprev;  
    double *u_end  = g->u_end;  
    double *u_f    = g->u_f;
    double cdist, fdist;

    if(debug && (myid==0)){
        res_loc = sqrt(res_loc);
        vector_dist(g, u_c,   u_cprev,&cdist); 
        vector_dist(g, u_end, u_f,    &fdist); 
        printf("kpar:%d, myid:%d, cvdist:%13.8e, fvdist:%13.8e , loc_res:%13.8e\n", kpar, myid, cdist, fdist, res_loc); 
    }
    // Fine solver has already evolved over the time slice, local residual should be ZERO
    if(mytid == kpar-1){ 
        if(res_loc > 0 ) {
            res_loc = sqrt(res_loc);
            vector_dist(g, u_c,   u_cprev, &cdist); 
            vector_dist(g, u_end, u_f,     &fdist); 
            WARN("Local residual should be ZERO, but it's NOT, details : \n \t kpar:%d, myid:%d, cvdist:%13.8e, fvdist:%13.8e, loc_res:%13.8e . \n\n", kpar, myid, cdist, fdist, res_loc); 
        }
    }

    if(myid==numprocs-1){
        if(pipelined == 0)
	        printf("kpar:%d, myid:%d, max_res:%13.8e \n", kpar, myid, max_res);
        else 
	        printf("kpar:%d, myid:%d, local_res:%13.8e \n", kpar, myid, res_loc);
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
   
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if(myid==0){  
        va_list args;
        va_start(args, fmt);
        fputs("ERROR: ", stderr);
        vfprintf(stderr, fmt, args);
        va_end(args);
    }
    MPI_Abort(MPI_COMM_WORLD, 2018); // the project since 2018 new year

    exit(1);
}

