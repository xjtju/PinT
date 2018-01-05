#include "Driver.h"
#include "HeatGrid.h"
#include "HeatSolverF.h"

/**
 * NOTE : the ini configuration file as the first command line parameter
 */
int main(int argc, char* argv[]) {
    
    //load global configuration and init MPI  
    Driver driver;
    driver.init(argc, argv);

    // get configuration information
    PinT* conf = PinT::instance();

    // create the grid/mesh and solver 
    Grid *g = new HeatGrid(conf);
    g->init();

    g->guardcell(g->u_end);

    //g->output();
    //driver.Abort("高次元テスト2D HEAT:%d\n", 5);

    Solver *F = new HeatSolverF(conf,g);    
    Solver *G = new HeatSolverC(conf,g);

    // run the parareal algorithm 
    driver.evolve(g, G, F);

    // output result
    g->output_local(g->u_end, false);
    //g->output_global();
    
    driver.finalize();  

    delete F;
    delete G;
    delete g;

    return 0;
}

