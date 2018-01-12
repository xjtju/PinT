#ifndef PinT_OUTPUT_H
#define PinT_OUTPUT_H 1 

#ifdef _HDF5_
#include "hdf5.h"
#endif

#include <stdio.h>
#include "Grid.h"

// the Output is eventually indepent out of Grid since the 2018 new year  
// because the output task become more and more heavier.
//
// WARN : 
//   the class is not well designed in the current version, only used for debug when the data is very small.
//   HDF5 support is added for large data,  but only for each local grid.
//
// NOTE : the caller must be responsible to judge whether the data is with or without border 
// and choose the proper output function.   

class Output {
private:
    Grid *grid;
    int ndim; 
    int sx, sy, sz; 
    int nx, ny, nz;
    int nguard;
    int idx, idy, idz;

    int idz_; // other grid's idz
    int idy_;
    int idx_;
    int spnumz;

public:

    Output(Grid *g){
        this->grid = g;
        ndim = g->ndim;
        sx = g->sx;
        sy = g->sy;
        sz = g->sz;
        
        nx = g->nx;
        ny = g->ny;
        nz = g->nz;

        idx = g->idx;
        idy = g->idy;
        idz = g->idz;
        idx_ = idx;
        idy_ = idy;
        idz_ = idz;

        nguard = g->nguard;

        spnumz = g->spnumz;
    }

    inline void coord(FILE *fp, int *cds) {
        fprintf(fp, "grid coords : [ %d ", cds[0]); 
        idx_=cds[0]*nx;
        if(ndim>=2){
            fprintf(fp, ", %d ", cds[1]); 
            idy_=cds[1]*ny;
        }
        if(ndim>=3) { 
            fprintf(fp, ", %d ", cds[2]);
            idz_=cds[2]*nz;
        } 
        fprintf(fp, "]\n"); 
    }

    // result output and debug 
    // NOTE : 
    // because the framework is mainly used for PinT performance testing,
    // the result is not important at current stage, so formal result output function is not yet provided,
    // such as HDF5 format output etc.  
   
    // in the current version, the global idz is not well supported, 
    // in order to avoid confusing, when outputing the global data, it is better to set withIdz=false
    // if set the withIdz=true, make sure to call output.coord() first to set the global z index correctly
    void var_inner_Z(FILE *fp, double *p, bool withIdz);  

    // output outer mainly used for debug
    void var_outer(FILE *fp, double *p); 
    void var_outer_X(FILE *fp, double *p); 
    void var_outer_Y(FILE *fp, double *p); 
    void var_outer_Z(FILE *fp, double *p); 

    //write all the solution data to HDF5
    int write_h5(char *fname, double *p);
};

#endif
