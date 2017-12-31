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
    g->create_topology();
    
    driver.Abort("高次元テスト:%d\n", 2);

    Solver *F = new HeatSolverF(conf,g);    
    Solver *G = new HeatSolverC(conf,g);
   
    // run the parareal algorithm 
    driver.evolve(g, G, F);

    // output result
    if(driver.myid == driver.numprocs-1) { 
        for(int i=0; i<g->size; i++){
            printf("%f\n", g->u_end[i]);
        }
        printf("kpar=%d\n",driver.kpar);
    }
    
    driver.finalize();  

    delete F;
    delete G;
    delete g;

    return 0;
}

