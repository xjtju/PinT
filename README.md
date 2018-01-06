# PinT Performance Testing Framework

A performance and convergency testing framework for Parallel-in-Time methods, especially for [Parareal](https://en.wikipedia.org/wiki/Parareal).

## Description

The framework is heavily inspired by Seigo Imamura's master dissertation "Hybrid-Parallel Performance Evaluation of Domain Decomposition and Parareal methods for Diffusion Problem" at Kobe University, Japan. It aims at facilitating the performance test of time-space hybrid parallel methodology in massive scientific computing or numerical simulation domain, especially Parareal,which is obtained more and more attention and research in the recent ten years since it was fully introduced by Martin J. Gander and Stefan Vandewalle in their paper "analysis of the parareal time-parallel time-integration method" at 2007.

The framework has implemented the Parareal algorithm skeleton in an uniform mesh, and also a linear system solver based on biconjugate gradient stabilized method, for running a test, the users only need to provide problem-specific stencil code, for example, the 2D heat equation's 5-point stencil. All the parameters controlling the time-space domain division, convergency, time-step etc. can be  predefined through an .INI file, easily be changed and tuned. If the default funtionality of it cannot be able to support some specific problem, it can be easily be extended by writting a new implementation in problem-specific sub classes. 

The framework is mainly written by C++ for good template and extension, most BLAS related calculations is performed by Fortran for performance reason and easy matrix manipulation. It is very light-weight, the only third library it used is [inih](https://github.com/benhoyt/inih), a small but excellent .INI file parser. 

At the current stage, only 1D or 2D can be automatically supported, 3D functions has being developed, maybe after two weeks. we plan to use [PMlib](https://github.com/avr-aics-riken/PMlib) as the performance moniter.  

## Design and implementation 

There are four main objects in the framework.
- PinT: for ini configuration, see the pint.ini sample
- Grid: managing the unifor mesh, holding the physical variables, iniializing variables and automaticaly applying boundary conditions, synchonizing guard cells, output result, and so on. 
- Driver: the driver of the PinT process. it implements the Parareal algorithm, and controls the execution flow of the program, especially the time parallel flow. The core function of the Driver is evolve(), it provides the template of Parareal algorithm, and drives the problem-specific coarse/fine solvers to evolve along the time slices within the whole time domain.    
- Solver: the abstract interface of all solvers for the framework, its most imporant sub class is linear solver based on [PBiCGStab](https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method).

## DEMO

For example, 2D heat equation, sample codes are under the src/heat directory.

1. establishs a grid (HeatGrid) for heat diffuse problem, the main and only task is to set the iniialization value to variables.
2. creates a sub class (HeatSolver) of the default linear solver (PBiCGStab), performs the stencil operations, the sample code uses the classic Crank–Nicolson method for deducing the stencil of heat equation.
3. defines fine/coarse solver based on the HeatSolver, and the fine and coarse solver is not necessary to use the same linear solver and time integrating method. 
4. combines the HeatGrid and the fine/coarse solver in the main program.
5. changes the .ini file according to the real run time envirement and the test request.

The main program sample for heat equation:

```c++
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

    // output result to disk for debug or post-processing 
    g->output_local(g->u_end, false); // for each process, create a unique file
    g->output_global();  // aggregates the result from all the process, and dumps out one file.
    
    driver.finalize();  // quit MPI 

    delete F;  
    delete G;
    delete g;  //free memory

    return 0;
}

```

