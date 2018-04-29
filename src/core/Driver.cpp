#include "Driver.h"
#include "blas.h"

#define TRACE(fmt, ...) ( (conf->verbose && mysid==0) ? printf(fmt, __VA_ARGS__ ): 0 ) 

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

/* 
 * the PinT algorithm template
 *
 * SKIP MODE:
 *   IF fine solver has already continuously evolved over the time slice, denoted by TS_{k}, 
 *   it is no meaningful to repeat the same calcaluatons until the TS_{k}, 
 *   but all the time slices before TS_{k} have to wait for convergency check, 
 *   there is NO way to quit safely for saving computing resources according to current MPI specification.
 *   The only thing we can do is to skip these useless send/receive operations to reduce a few of communication overhead.  
 *   Of course, the fine/coarse solver's calculation can also be skipped safely, but the MPI process has to wait anyway.
 *   At current, we only skip send/receive, the calculation skip is tested, but commented out. 
 * 
 * NOTE:
 *   the guardcell synchonization is the task of space parallel (Grid), not time parallel (Driver) 
 *   according to the design principle of our framework. 
 *   Guardcell synchonization is not necessary to be called explictly here except the end of the time loop    
 *   for the consistency of output result with guardcell.
 *
 *   There should be several different forms of the residual expressions,
 *   the current version used a relative value |res|^2 / |u|^2 
 *   it is recommended to choose a proper calcaluation formula for specific problem.   
 */
void Driver::evolve(Grid* g, Solver* G, Solver* F){
    // output the initial values for debugging or ...
    if(conf->dump_init) {
        g->output_global_h5("ini");
        INFO("\ninitial data has been dumped to ...all.ini.h5 \n\n");
    }
    // for convience only, set the pointer to grid inner variables
    double *u_start= g->u_start; 
    double *u_end  = g->u_end;  
    double *u_c    = g->u_c;  
    double *u_f    = g->u_f; 

    double *sendslns = G->sendslns();
    double *recvslns = G->recvslns();
    size_t solnsize  = G->solnsize();

    double relax_factor = conf->relax_factor;  // convergence accelerator 

    int source, dest, tag;
    size_t size = g->size; 
    int ierr;
    MPI_Status  stat;
    unsigned long icount = 0;  

    // except the first time slice, all others need to receive U^{0}_{n-1} as its start value  
    monitor.start(Monitor::RECV);
    if(isRecvSlice(1)) {
        source = myid - spnum;
        tag = myid*100;
        MPI_Recv(recvslns, solnsize, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &stat);
        G->unpack();
    } 
    monitor.stop(Monitor::RECV);

    //coarse 
    monitor.start(Monitor::CSolver);
    blas_cp_(u_c, u_start, &size);  
    icount = G->evolve();
    // unify the process flow of send/recv, in a sence, u_c can be regarded as u_end at the init step 
    blas_cp_(u_end, u_c, &size);  
    monitor.stop(Monitor::CSolver, 1, icount);
    
    monitor.start(Monitor::SEND);
    // except the last time slice, all others need to send the coarse(estimate) value to its next slice  
    if(isSendSlice(1)){
        dest = myid + spnum;
        tag  = (myid + spnum)*100;
        //send to next time slice
        G->pack(); //pack is necessary only at initialization of parareal
        ierr = MPI_Rsend(sendslns, solnsize, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);    
    }
    monitor.stop(Monitor::SEND);

    double res_loc, res_sp, max_res, nrm_loc, nrm_sp;
    kpar = 0;
    for(int k=1; k<=conf->kpar_limit; k++)
    {
        res_loc = 0.0;
        res_sp  = 0.0;
        max_res = 0.0;
        nrm_loc = 0.0;
        nrm_sp  = 0.0;

        kpar = kpar + 1;
        
        // step1: backup previous solutions of coarse solver
        //blas_cp_(u_cprev, u_c, &size); //this step is not necessary at the following of fine solver 
        G->backup_prevs();

        // step2: fine solver parallel run based on U^{k-1}_{n-1}
        monitor.start(Monitor::FSolver);
        blas_cp_(u_f, u_start, &size);
        //if(!isSkip(k))
        icount = F->evolve(); 
        //else TRACE("%d, F is skiped\n",mytid);
        monitor.stop(Monitor::FSolver, 1, icount);

        // step3: receive the latest solutions from the previous slice
         monitor.start(Monitor::RECV); 
        if(isRecvSlice(k)){
	        source = myid - spnum;
	        tag    = myid*100 + kpar;
	        MPI_Recv(recvslns, solnsize, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &stat);
        } else TRACE("SKIP RECV : t=%d/s=%d \n", mytid, mysid); 
        if(mytid!=0) G->unpack();
        monitor.stop(Monitor::RECV);

        // step4: run coarse solver based on the latest solutions just received 
        monitor.start(Monitor::CSolver);
        blas_cp_(u_c, u_start, &size); 
        //if(!isSkip(k)) 
        icount = G->evolve();
        //else TRACE("%d, G is skiped\n",mytid);
        monitor.stop(Monitor::CSolver, 1, icount);

        // step5: correct the value and calculate the local residual
        //if(!isSkip(k)) {
        pint_sum(g, &conf->num_std, sendslns, F->curr_solns(), G->curr_solns(), G->prev_solns(), &relax_factor, &res_loc, &nrm_loc);  
        G->update_uend(); //make sure the u_end has the latest solution 
        //}else TRACE("%d, SUM is skiped\n",mytid);
        
        // step6: send to latest correced solutions to the next slice
        monitor.start(Monitor::SEND); 
        if(isSendSlice(k)){
	        dest = myid + spnum;
	        tag  = (myid + spnum)*100 + kpar;
            //here, pack is unnecessary, pint_sum() do it
	        ierr = MPI_Send(sendslns, solnsize, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);    
        } else TRACE("SKIP SEND : t=%d/s=%d \n", mytid, mysid); 
        monitor.stop(Monitor::SEND); 

        //STEP7 : gather residual
        monitor.start(Monitor::RES);
        g->sp_allreduce(&res_loc, &res_sp);
        g->sp_allreduce(&nrm_loc, &nrm_sp);
        //printf("res_loc=%e, res_sp=%e\n", res_loc, res_sp);
        if( nrm_sp < 1.0e-108 ) nrm_sp = smlr; // avoid divided by ZERO
        res_sp = res_sp/nrm_sp;
        if( res_sp < smlr )  res_sp = 0.0;
        //MPI_Allreduce(&res_loc, &res_sp,  1, MPI_DOUBLE, MPI_SUM, sp_comm);
        if(pipelined == 0) {
            g->allreduce(&res_sp, &max_res, MPI_MAX);  // time space is enough ?!
            //MPI_Allreduce(&res_sp,  &max_res, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        }
        max_res = sqrt(max_res);
        monitorResidual(g, res_loc, max_res, size);
        monitor.stop(Monitor::RES);

         //STEP8 : check the convience
        if(pipelined == 0){
            if( max_res < conf->converge_eps) {
                if(myid==numprocs-1) printf("Parareal is converged at %d iteration, res = %e .\n", k , max_res);
                break;   
            }
        }
    }
    g->guardcell(u_end);
    MPI_Barrier(MPI_COMM_WORLD);
}

void Driver::finalize(bool pfile) {
    monitor.gather();

    if(pfile) {
        char fname[40];
        memset(fname, 0, sizeof(char)*40);
        sprintf(fname, "%s_%s_basic.txt", conf->monitor_pre, jobid); 
        monitor.print(fname, "The TAO of Programming", "The PinT performance test framework");
    
        memset(fname, 0, sizeof(char)*40);
        sprintf(fname, "%s_%s_detail.txt", conf->monitor_pre, jobid); 
        monitor.printDetail(fname);
    }

    //Sometimes, profiling information cannot be completely written into common files like above in HPC environments,
    //but stdout/stderr has no problem.  
    monitor.print(stderr, "ITO supercomputer in Kyushu University", "The PinT performance test framework");
    monitor.printDetail(stderr);

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
    double cdist;

    if(debug && (mysid==0)){
        res_loc = sqrt(res_loc);
        vector_dist(g, u_c,   u_cprev,&cdist); 
        printf("kpar:%d, mytid:%d, cvdist:%13.8e, loc_res:%13.8e\n", kpar, mytid, cdist, res_loc); 
    }
    // Fine solver has already evolved over the time slice, local residual should be ZERO
    // mytid starts from ZERO, kpar starts from ONE
    if(mytid == kpar-2){ 
        if(res_loc > 0 ) {
            res_loc = sqrt(res_loc);
            vector_dist(g, u_c,   u_cprev, &cdist); 
            WARN("Local residual should be ZERO, but it's NOT, details : \n \t kpar:%d, mytid:%d, myid:%d, cvdist:%13.8e, loc_res:%13.8e . \n\n", kpar, mytid, myid, cdist, res_loc); 
        }
    }

    if(myid==numprocs-1){
        if(pipelined == 0)
	        printf("kpar:%03d, mytid:%03d, max_res:%13.8e \n", kpar, mytid, max_res);
        else 
	        printf("kpar:%03d, mytid:%03d, local_res:%13.8e \n", kpar, mytid, res_loc);
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

    //if(myid != 0) return ;
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

    va_list args;
    va_start(args, fmt);
    fputs("ERROR: ", stderr);
    vfprintf(stderr, fmt, args);
    va_end(args);
    
    MPI_Abort(MPI_COMM_WORLD, 2018); // the project since 2018 new year

    exit(1);
}

