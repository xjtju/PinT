#ifndef PinT_GRID_H
#define PinT_GRID_H 1 

#include "common.h"

/**
 * in current, the class is only for holding the physical variables.
 * In the further, it can be extended to a mesh structure
 */

class Grid {

public:
    //the mesh size of the whole space grid
    int  nx;
    //inline void set_nx(int n){
    //    this->nx = n;
    //}

    int nguard = 1;
    //inline void set_nguard(int n){
     //   this->nguard = n;
    //}

    int size;
    double dx,dt;   

    // the physical variables     
    double *x;

    Grid(int nx, int ng, double dx, double dt); 
    ~Grid();

    virtual int init()=0;

    //boundary condition
    virtual void bc()=0;
};
#endif
