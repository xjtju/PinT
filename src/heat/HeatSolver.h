#ifndef _PinT_HEATSOLVER_H_
#define _PinT_HEATSOLVER_H_ 1

#include "heat.h"
#include "Solver.h"
#include "PBiCGStab.h"

/**
 * 1D/2D/3D heat equation with Crank-Nicolson
 *
 * The class is also a template for linear system which can be directly integrated in time loops
 */
class HeatSolver : public Solver {
private:
    void setup();


public:
    
    double k ;  // diffuse coefficient

    // lamda = k*dt/dx^2   
    // when dx, dy, dz are not equal with each other, the lamdas are unequall neither.
    double lamda_x;
    double lamda_y;
    double lamda_z;
    double &lamda = lamda_x;

    double lamdaxyz[3];
    
    double *soln; // the pointer to grid->u_f/u_c

    HeatSolver(PinT *c, Grid *g); 
    HeatSolver(PinT *c, Grid *g, bool isFS);


    virtual ~HeatSolver() {
        delete hypre;
    }
    
    void evolve();         // evolve over a time slice 

    inline void stencil() {
        if(ndim==3)      stencil_heat_3d_(grid->nxyz, lamdaxyz, &nguard, soln, A);
        else if(ndim==2) stencil_heat_2d_(grid->nxyz, lamdaxyz, &nguard, soln, A);
        else if(ndim==1) stencil_heat_1d_(grid->nxyz, lamdaxyz, &nguard, soln, A);
    }

    inline void rhs() {
        if(ndim==3)      rhs_heat_3d_(grid->nxyz, lamdaxyz, &nguard, soln, b);
        else if(ndim==2) rhs_heat_2d_(grid->nxyz, lamdaxyz, &nguard, soln, b);
        else if(ndim==1) rhs_heat_1d_(grid->nxyz, lamdaxyz, &nguard, soln, b);
    }

};

#endif
