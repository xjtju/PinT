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
 */
#ifdef _TEST_HEAT_

#include "Driver.h"
#include "HeatGrid.h"
#include "HeatSolverF.h"

int main(int argc, char* argv[]) {
    
    //load global configuration and init MPI  
    Driver driver;
    driver.init(argc, argv);

    // get init-ready system information
    PinT* conf = PinT::instance();

    // create the grid/mesh and solver 
    Grid *g = new HeatGrid(conf);
    g->init();


    Solver *F = new HeatSolverF(conf,g);   // fine solver 
    Solver *G = new HeatSolverC(conf,g);   // coarse solver

    // run the parareal algorithm 
    driver.evolve(g, G, F);

    // output result to disk and for post-processing 
    g->output_local(g->u_end, true);
    g->output_global_h5("heat");
    
    //driver.Abort("高次元テスト3D HEAT:%d\n", 3); // DEBUG
    driver.finalize();  // quit MPI 

    delete F;  
    delete G;
    delete g;  //free memory

    return 0;
}
#endif
