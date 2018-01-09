#ifndef PinT_OUTPUT_H
#define PinT_OUTPUT_H 1 

#include <stdio.h>
#include "Grid.h"

// the Output is indepent out of Grid from the 2018 new year  
// because the output task become more heavy 
class Output {
private:
    Grid *grid;
    
    int sx, sy, sz;
    int nguard;
    int idx, idy, idz;

public:

    Output(Grid *g){
        this->grid = g;

        sx = g->sx;
        sy = g->sy;
        sz = g->sz;

        idx = g->idx;
        idy = g->idy;
        idz = g->idz;

        nguard = g->nguard;
    }

    void var_outer(FILE *fp, double *p); 
    void var_outer_X(FILE *fp, double *p); 
    void var_outer_Y(FILE *fp, double *p); 
    void var_outer_Z(FILE *fp, double *p); 
};

#endif
