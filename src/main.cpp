#include "Driver.h"
#include "HeatGrid.h"
#include "HeatSolverF.h"

/**
 * the example program for employing parareal method (PinT) to solve real problem based on the framework.
 *
 * NOTE : the ini configuration file as the first command line parameter, 
 * default configuration file is "pint.ini", under the same directory with the executable file. 
 */
int main(int argc, char* argv[]) {
    
    //load global configuration and init MPI  
    Driver driver;
    driver.init(argc, argv);

    // get init-ready system information
    PinT* conf = PinT::instance();

    // create the grid/mesh and solver 
    Grid *g = new HeatGrid(conf);
    g->init();
    g->guardcell(g->u_end);
    g->bc(g->u_end);
    g->output_local(g->u_end, false);

    driver.Abort("高次元テスト3D HEAT:%d\n", 1); // DEBUG

    Solver *F = new HeatSolverF(conf,g);   // fine solver 
    Solver *G = new HeatSolverC(conf,g);   // coarse solver

    // run the parareal algorithm 
    driver.evolve(g, G, F);

    // output result to disk and for post-processing 
    g->output_local(g->u_end, false);
    g->output_global();
    
    driver.finalize();  // quit MPI 

    delete F;  
    delete G;
    delete g;  //free memory

    return 0;
}

