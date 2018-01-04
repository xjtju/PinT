#ifndef _PinT_HEATSOLVER_F_H_
#define _PinT_HEATSOLVER_F_H_

#include "HeatSolver.h"

class HeatSolverF : public HeatSolver {
public:
    HeatSolverF(PinT *conf, Grid *g);
    
    double* fetch();
};

class HeatSolverC : public HeatSolver {
public:
    HeatSolverC(PinT *conf, Grid *g);
    double* fetch();
};

#endif
