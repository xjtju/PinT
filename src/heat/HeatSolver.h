#ifndef _PinT_HEATSOLVER_H_
#define _PinT_HEATSOLVER_H_ 1

#include "heat.h"
#include "DefaultSolver.h"

/**
 * 1D/2D/3D heat equation with Crank-Nicolson
 *
 * The class is an example about how to use DefaultSolver as the time integration template
 */
class HeatSolver : public DefaultSolver {
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
    
    HeatSolver(PinT *c, Grid *g); 
    HeatSolver(PinT *c, Grid *g, bool isFS);


    virtual ~HeatSolver() { } 

    void stencil() {
        if(ndim==3)      stencil_heat_3d_(grid->nxyz, lamdaxyz, &nguard, soln, A);
        else if(ndim==2) stencil_heat_2d_(grid->nxyz, lamdaxyz, &nguard, soln, A);
        else if(ndim==1) stencil_heat_1d_(grid->nxyz, lamdaxyz, &nguard, soln, A);
    }

    void rhs() {
        if(ndim==3)      rhs_heat_3d_(grid->nxyz, lamdaxyz, &nguard, soln, b);
        else if(ndim==2) rhs_heat_2d_(grid->nxyz, lamdaxyz, &nguard, soln, b);
        else if(ndim==1) rhs_heat_1d_(grid->nxyz, lamdaxyz, &nguard, soln, b);
    }

};

#endif
