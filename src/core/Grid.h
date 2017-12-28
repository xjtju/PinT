#ifndef PinT_GRID_H
#define PinT_GRID_H 1 

#include "common.h"
#include "PinT.h"
/**
 * in current, the class is only for holding the physical variables used by PinT.
 * In the further, it can be extended to a mesh structure
 *
 * The iniialization of grid variables and boundary condition is done by Grid, so do the guard cell synchonization.
 *   
 *  NOTE :
 *    In PinT, for each physical variable, it has to use five times memory space as its real size. 
 *
 *    There should be only one instance of Grid in the whole application ,that is singleton.
 *    
 */

class Grid {

public:
    PinT *conf;

    //the mesh size of the whole space grid
    int  nx;

    int nguard = 1;

    int size;
    double dx;   

    
    // the output of fine/coarse solver      
    double *u_f;  
    double *u_c;
    double *u_cprev; // the structure holder for the coarse solver at the previous iteration of the the same slice 
    
    // the physical variables or the solution at the edge of the time slice   
    double *u_start; // the latest solution of the current time slice start point or the previous slice end 
    double *u_end;   // the latest solution of the current time slice end point or the next slice start 
    double *u;       // pointing the same variable with u_end


    Grid(PinT* conf); 
    ~Grid();

    virtual int init()=0;

    //boundary condition
    virtual void bc()=0;

    void guardcell() {};
};
#endif
