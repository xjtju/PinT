#ifndef _PinT_HEATSOLVER_F_H_
#define _PinT_HEATSOLVER_F_H_

#include "HeatSolver.h"

class HeatSolverF : public HeatSolver {
public:
    HeatSolverF(PinT *conf, Grid *g):HeatSolver(conf,g){

        this->steps = conf->f_steps;  
    }

    double* fetch(){
        return grid->u_f;
    }
};

class HeatSolverC : public HeatSolver {
public:
    HeatSolverC(PinT *conf, Grid *g):HeatSolver(conf,g){

        this->steps = conf->c_steps;  
    }
    double* fetch(){
        return grid->u_c;
    }
};

#endif
