# PinT Performance Testing Framework

A performance and convergency testing framework for Parallel-in-Time methods, especially for [Parareal](https://en.wikipedia.org/wiki/Parareal).

## Description

The framework is heavily inspired by Seigo Imamura's master dissertation *"Hybrid-Parallel Performance Evaluation of Domain Decomposition and Parareal methods for Diffusion Problem"* at Kobe University, Japan. It aims at facilitating the performance test of time-space hybrid parallel methodology in massive scientific computing or numerical simulation domain, especially Parareal,which is obtained more and more attention and research in the recent ten years since it was fully introduced by Martin J. Gander and Stefan Vandewalle in their paper *"analysis of the parareal time-parallel time-integration method"* at 2007.

The framework has implemented the **Parareal algorithm skeleton** in an uniform mesh, and also a linear system solver based on biconjugate gradient stabilized method, for running a test, the users only need to provide problem-specific stencil code, for example, the 2D heat equation's 5-point stencil. All the parameters controlling the time-space domain division, convergency, time-step etc. can be  predefined through an .INI file, easily be changed and tuned. If the default funtionality of it cannot be able to support some specific problem, it can be easily be extended by writting a new implementation in problem-specific sub classes. 

The framework is mainly written by C++ for good template and extension, most BLAS related calculations is performed by Fortran for performance reason and easy matrix manipulation. It is very light-weight, the only third library it used is [inih](https://github.com/benhoyt/inih), a small but excellent .INI file parser. 

At the current stage, only 1D or 2D can be automatically supported, 3D functions has being developed, maybe after two weeks. we plan to use [PMlib](https://github.com/avr-aics-riken/PMlib) as the performance moniter.  

## Design and implementation 

There are four main objects in the framework.
- **PinT**: for ini configuration, see the pint.ini sample file for details
- **Grid**: manages the uniform mesh, holds the physical variables, iniializes variables and applys boundary conditions, synchonizes guard cells, outputs result, and so on. 
- **Driver**: the driver of the PinT process. it implements the Parareal algorithm, and controls the execution flow of the program, especially the time parallel flow. The core function of the Driver is evolve(), it provides the template of Parareal algorithm, and drives the problem-specific coarse/fine solvers to evolve along the time slices within the whole time domain.    
- **Solver**: the abstract interface of all solvers for the framework, its most imporant sub class is linear solver based on [PBiCGStab](https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method).

## DEMO

For example, [1D/2D heat equation](https://commons.wikimedia.org/wiki/File:Heatequation_exampleB.gif), sample codes are under the src/heat directory.

1. establishs a grid (HeatGrid) for heat diffuse problem, the main and only task is to set the initial value to variables.

```c++
HeatGrid::HeatGrid(PinT *conf) : Grid(conf){ } 

void HeatGrid::init(){
    long ind = 0;
    double x, unk;
    for(int i = nguard; i<nx+nguard ; i++){
        ind = i;
        x = this->getX();    // global coordincate of the cell
        unk = cos(2*x);      // set the initial temperature
        // set the variables used by Parareal method 
        u_start[ind] = unk;  // for start point of the current time slice
        u_f[ind] = unk;      // for fine solver, not necessary, it will be also set automatically   
        u_c[ind] = unk;      // for coarse solver, not necessary, it will be also set automatically    
    }
}
```
2. creates a sub class (HeatSolver) of the default linear solver (PBiCGStab), performs the stencil operations, the sample code uses the classic [Crank–Nicolson method](https://en.wikipedia.org/wiki/Crank%E2%80%93Nicolson_method) for deducing the stencil of heat equation.

```c++
// set diffuse coefficient and tune the default parameter, problem specific
void HeatSolver::setup(){
    this->eps = 1.0e-6;  // change the default value of the super class
    k = 0.061644;        // diffuse coefficient 
}

// 1D, not used Fortran

//calcaluate the residual r = b - Ax
void HeatSolver::cg_rk1d(double *r, double *x, double *b){
    for(int i=nguard; i<nx+nguard; i++){
        double ax = -lamda*x[i-1] + (1+2*lamda)*x[i] - lamda*x[i+1];  
        r[i] = b[i] - ax;
    }
}

// matrix * vector,  the stencil matrix , v = Ay
void HeatSolver::cg_Xv1d(double* v, double *y) {
    for(int i=nguard; i<nx+nguard; i++){
        v[i] = -lamda*y[i-1] + (1+2*lamda)*y[i] - lamda*y[i+1];
    }
}

// calcaluate b of Ax=b,  RHS
void HeatSolver::cg_b1d(double *x){
    for(int i=nguard; i<nx+nguard; i++){
        b[i] = lamda*x[i-1] + (1-2*lamda)*x[i] + lamda*x[i+1]; 
    }
}

```
3. defines fine/coarse solver based on the HeatSolver, and the fine and coarse solver is not necessary to use the same linear solver and time integrating method. For the Fine and Coarse solver, the only thing is to set their specific variables in most cases. 

```c++
HeatSolverF::HeatSolverF(PinT *conf, Grid *g):HeatSolver(conf,g, true) {
    lamda = k * conf->f_dt / (2*g->dx*g->dx); // provides the lamda's value
}

HeatSolverC::HeatSolverC(PinT *conf, Grid *g):HeatSolver(conf,g, false){
    lamda = k * conf->c_dt / (2*g->dx*g->dx);
}
```

4. combines the HeatGrid and the fine/coarse solver in the main program.

```c++
int main(int argc, char* argv[]) {
    
    //load global configuration and init MPI  
    Driver driver;
    driver.init(argc, argv);

    // get init-ready system information
    PinT* conf = PinT::instance();

    // create and init the grid/mesh and solver 
    Grid *g = new HeatGrid(conf);
    g->init();
    Solver *F = new HeatSolverF(conf,g);   // fine solver 
    Solver *G = new HeatSolverC(conf,g);   // coarse solver

    // run the parareal algorithm 
    driver.evolve(g, G, F);

    // output result to disk for debug or post-processing, the two steps is not necessary for performance test. 
    g->output_local(g->u_end, false); // for each process, it will create a unique file
    g->output_global();  // aggregates the result from all the process, and dumps out to one file.
    
    driver.finalize();  // quit MPI 

    delete F;  
    delete G;
    delete g;  //free memory

    return 0;
}

```
5. changes the .ini file according to the real run time envirement and the test request. See pint.ini sample for details, and the .INI file is very direct and simple. 



From the above sample codes, in most cases it is not necessary for users to care many boilerplate tasks explicitly such managing MPI envirement, mesh division, guard cell synchonization etc. The framework can perform most housekeeping tasks well, so users can focus on phyical moddel or problem itself.  

## Notice 

The framework is designed for performance test originally, not for the high accuracy scientific computing or numerical simulations. There are not any assumptions about the measurement unit standard for all the physical variables used in the program. You must be careful to make sure the physical UNIT consistency for specific problem.  

In the current version, in order to make the framework usable as earlier as possible to verify our design ideas, no special get/set methods are provided for most grid variables and configuration parameters. You must be careful when directly changing these values in some cases.     

## My consideration about Parareal

1. It uses more system memory (perhaps 5 times) than the traditional pure space parallel methodology. 
2. In practice, for many real-world multi-physical simulations, it is hard to pack all the related solvers into Parareal algorithm framework. How to do ?   

