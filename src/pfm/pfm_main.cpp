#ifdef _TEST_PFM_

#include "Monitor.h"
#include "Driver.h"

#include "CNSolver.h"
#include "BD4Solver.h"
#include "EulerSolver.h"
#include "PFMGrid.h"
#include "PFMParams.h"

/**
 * the example program for employing parareal method (PinT) to solve real problem based on the framework.
 *
 * NOTE : 
 *
 * (1). the ini configuration file is as the 1st command line parameter, 
 * default configuration file is "pint.ini", under the same directory with the executable file.
 *
 * (2). if using profiling, the profiling output file is as the 2nd command line parameter, 
 * default profiling file is [monitor_pre].000.000.txt
 *
 * OOPS: 
 *   I have encountered the 'super-linear speedup' phenomenon at the ITO supercomputer of Kyusyu university.
 *   When the  grid size is fixed at 128^3, the performance of 256 cores is nearly 3-fold of 128 cores.  
 *   The CPU is Xeon Gold 6154(Skylake-SP) : 18 cores per CPU, non-inclusive 1.375 MB LLC per core. 
 *  
 *   The memory for holding solution variables ~ 128^3*8 bytes ~ 16MB.   
 *   There are several data structures used by parareal algorithm and linear solvers,
 *   if both fine and coarse solver using the combination of CN and PBiCG,  approximately 16MB*23, 
 *   for each core is 16MB*23/256 ~ 1.4375MB, just a bit more than 1.375MB.
 *
 *   Perhaps I think it can be as an explanation to the best of my knowledge.  
 */

int main(int argc, char* argv[]) {
    
    //load global configuration (parameters) and init MPI  
    Driver driver;
    driver.init(argc, argv);

    // get init-ready system information
    PinT* conf = PinT::instance();
    
    // create the grid/mesh 
    Grid *g = new PFMGrid(conf);
    g->init();
    
    // choose and setup the proper coarse/fine solver
    PFMParams param;
    param.init(conf,g);
    Solver *G, *F;
    
    if(param.csolver == param.ID_CN) 
         G = new CNSolver (conf, g, false);    
    else G = new BD4Solver(conf, g, false);  

    if(param.fsolver == param.ID_EU) 
         F = new EulerSolver(conf, g, true);
    else F = new CNSolver   (conf, g, true);   

    //driver.Abort("高次元テスト3D PFM:%d\n", 5); // DEBUG
    // run the parareal algorithm 
    driver.evolve(g, G, F);

    // output result to disk and for post-processing 
    // For large scale performance test, it is best to close the output
    //g->output_local(g->u_end, true); // ASCII mode
    g->output_global_h5("pfm");        // HDF5  mode 
    
    driver.finalize();  // quit MPI 

    delete F;  
    delete G;
    delete g;  //free memory

    return 0;
}

#endif
