
**| Framework and tools** 

- include : header files of the framework
- .py .sh  : plot scripts and job running scripts

- core : the main template codes of the framework, for higher abstraction 
  - PinT   : parameter and configuration management 
  - Grid   : manages space division and communication 
  - Driver : manages time division and time integration
  - Output : writes the result or other variables to disk, ASCII for small data and HDF5 for large data is supported.
  - LS     : the abstraction of linear solver
    - SOR       : Black-Red successive over relaxation method 
    - PBiCGStab : Biconjugate gradient stabilized method with precondictioner
  - Solver : the time integrator used by time parallel
    - NewtonSolver : the template of Newton-Raphson method
  - blas.f90   : common matrix calculations  
  - blascg.f90 : BiCG-specific calculations

- utils
  - util1d/2d/3d.f90 : fortran source file for multi-dimension data manipulation 
  - Monitor : A wrapper of PMLib for performance measurement
  - ini.c   : .INI parser from https://github.com/benhoyt/inih  


**| Examples and test cases**  

- heat : heat equation with constant diffuse coefficient, 1D/2D/3D
  - HeatGrid   : intializes variables
  - HeatSolver : a common time integrator for both coarse and fine solver 
  - HeatSolverF: example for coarse solver and fine solver
  - heat.f90   : the calculations of RHS and b for linear sytem Ax=b 
  - heat.ini   : configuration file
  - heat_main.cpp : main program entry


- pfm : phase field model - Allen-Cahn equation, 1D/2D/3D 
  - PFMGrid    : set initial value and provides customized BC function 
  - CNSolver   : Newton-Raphson time integrator using Crank-Nicolson finite difference method 
  - BD4Solver  : Newton-Raphson time integrator using 4th-order backward differentiation formula 
  - PFMParams  : the holder of the problem-specific parameters 
  - EulerSolver: explicit Euler method, unstable! just for example 
  - pfm.90/pfm_bd4.f90 : problem-specific matrix calculations
  - pfm.ini    : configuration file
  - pfm_main.cpp : main program entry

