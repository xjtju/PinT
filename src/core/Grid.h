#ifndef PinT_GRID_H
#define PinT_GRID_H 1 

#include "common.h"
#include "FortFunc.h"
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

    int ndim;

    //the mesh size of the whole space grid
    int  nxyz[3];
    int  nx;  //grid size without guard cells
    int  ny;
    int  nz;
    int inner_size; // size not including guard cells

    int nguard = 1;
    int bc_type;
    double bc_val;

    // in order to automatically adapt to multi-dimension  
    int ngxyz[3];  

    long size;  // size including guard cells
    int sxyz[3];
    int sx; //grid size with guard cells
    int sy;
    int sz;

    int gcsx; //guard cell size 
    int gcsy;
    int gcsz;

    double dx; //grid cell size   
    double dy;
    double dz;

    long idx; //global coord idx of the left-back-bottom point  
    long idy;
    long idz;

    // the output of fine/coarse solver      
    double *u_f;  
    double *u_c;
    double *u_cprev; // the structure holder for the coarse solver at the previous iteration of the the same slice 
    
    // the physical variables or the solution at the edge of the time slice   
    double *u_start; // the latest solution of the current time slice start point or the previous slice end 
    double *u_end;   // the latest solution of the current time slice end point or the next slice start 
    double *u;       // pointing the same variable with u_end
   
    int myid;   // MPI process id in global 
    int mysid;  // process id in space domain

    int spnum ; // the number of space partitions, parallel processes along the space domain 
    int spnumx;
    int spnumy;
    int spnumz;
    MPI_Comm *sp_comm; //space within the same time slice
    int st_rank;       // topology rank 
    MPI_Comm st_comm;  // topology comm
    int* coords;      // topology coordinates 
    int* coords_;     // used for hold coordinates of other grids 
    int* dims;
    int* periods;

    // for guardcell send and receive in one direction 
    double *gcell_sendx; // left & right
    double *gcell_sendy; // front & back
    double *gcell_sendz; // top & bottom

    double *gcell_recvx;
    double *gcell_recvy;
    double *gcell_recvz;

    
    // time slice information
    int mytid;
    int tsnum; 

    Grid(PinT* conf); 
    ~Grid();

    virtual int init()=0;

    // boundary condition
    void bc();
    void bc(double *d); 
    void bc_1d(double *d); 
    void bc_2d(double *d); 

    // create the grid topology
    void create_topology(); 
    void create_topology_1d();
    void create_topology_2d();

    // guardcell exchange  
    void guardcell();
    void guardcell(double *d);
    void guardcell_1d(double *d);
    void guardcell_2d(double *d);

    // grid topology. left:X:right; front:Y:back; top:Z:bottom
    int left, right, front, back, top, bottom;
  
    /**
     * get the XYZ global geographical value from their local outer index,
     * often used in iniialization and position-related calculations.
     */
    // ind is the outer index
    inline double getX(int ind){
        int gind = (idx+ind-nguard);
        return gind*dx + dx/2; 
    }
    inline double getY(int ind){
        int gind = (idy+ind-nguard);
        return gind*dy + dy/2;
    }
    inline double getZ(int ind){
        int gind = (idz+ind-nguard);
        return gind*dz + dz/2;
    }

    // get the linear index of the 1-2-3D index, the set of functions are not frequently used.
    inline long getInnerIdx(int ix){
        return ix-nguard;
    }
    inline long getInnerIdx(int ix, int iy){
        return nx*(iy-nguard)+(ix-nguard);
    }
    inline long getInnerIdx(int ix, int iy, int iz){
        return ny*nx*(iz-nguard)+nx*(iy-nguard)+(ix-nguard);
    }

    inline long getOuterIdx(int ix){
        return ix;
    }
    inline long getOuterIdx(int ix, int iy){
        return sx*iy + ix;
    }
    inline long getOuterIdx(int ix, int iy, int iz){
        return sy*sx*iz + sx*iy + ix;
    }
    
    // get rid of guard cells and pack inner grid data into a buffer  
    inline void pack_data(double *p, double *buf){
        switch(ndim) {
          case 1: pack_1d_(nxyz, &nguard, p, buf); break;
          case 2: pack_2d_(nxyz, &nguard, p, buf); break;
          case 3: printf("3D is not finished\n");  break;
        }
    }

    // space domain, for only one datum with double type and SUM operation 
    void sp_allreduce(double *d);  //d is input and output
    void sp_allreduce(double *d, double *o); //d is input, o is output
    // time-space domain
    void allreduce(double *d, double *o, int op);

    
    // result output and debug 
    // NOTE : 
    // because the framework is mainly used for PinT performance testing,
    // the result is not important at current stage, so formal result output function is not yet provided,
    // such as HDF5 format output etc.  

    // the basic output function, write grid local variable for debug, 
    // the caller must be responsible to judge whether the data is with or without border 
    // and choose the proper output function.   
    void output_var_inner(FILE *fp, double *data); 
    void output_var_outer(FILE *fp, double *data); 

    //the wrapper of output_var_*, writing the local variables into debug file only for the last time slice 
    void output_local(double *data, bool inner_only);

    // aggregate all the final results from all the grids within the same space domain and output 
    void output_global();

    inline void printf_coord(FILE *fp, int *cds) {
        fprintf(fp, "[ %d ", cds[0]); 
        if(ndim>=2) fprintf(fp, ", %d ", cds[1]); 
        if(ndim>=3) fprintf(fp, ", %d ", cds[2]); 
        fprintf(fp, "]\n"); 
    }

};
#endif
