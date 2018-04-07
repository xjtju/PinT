#ifndef PinT_GRID_H
#define PinT_GRID_H 1 

#include "common.h"
#include "FortFunc.h"
#include "PinT.h"
#include "Monitor.h"

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
 *  TO DO:
 *    As the number of parameters of grid increases, it is necessary to seperate the its configuration code from PinT.    
 *    
 */

class Grid {
 private:
    Monitor monitor;
    inline void dummy(const char* s) { if(myid==0) printf("WARN: dummy %s function is used  by GRID !\n", s); }

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

    int bc_type;  // 0: Dirichlet ;  1: homogeneous Neuman; 2: customize
    double bc_val;

    int bc_type_xl = 0;      // left x boundary  
    int bc_type_xr = 0;      // left x boundary  
    double bc_val_xl = 0.0;  
    double bc_val_xr = 0.0;  

    int bc_type_yl = 0;      // left x boundary  
    int bc_type_yr = 0;      // left x boundary  
    double bc_val_yl = 0.0;  
    double bc_val_yr = 0.0;  

    int bc_type_zl =0;      // left  z boundary (bottom/down) 
    int bc_type_zr =0;      // right z boundary (upper) 
    double bc_val_zl = 0.0;  
    double bc_val_zr = 0.0;  

    // in order to automatically adapt to multi-dimension  
    int ngxyz[3];  

    long size;  // size including guard cells
    int sxyz[3];
    int sx;     // grid size with guard cells
    int sy;
    int sz;

    int gcsx;  // guard cells size, NOTE : it is different with nguard
    int gcsy;   
    int gcsz;

    double dx; // grid cell size, geographical size   
    double dy;
    double dz;

    long idx; //for recording the global coord of the grid, start point is the most left-back-bottom point  
    long idy;
    long idz;

    // the temporary solution variables of fine/coarse solver      
    // or the data interface between Grid and Fine/Coarase Solver
    // before fine/coarse solver starts a new iteration, the Driver will copy the solution at starting point to u_f/u_c
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
    int* dims;        // topology's dim
    int* periods;     // topology's period

    // for guardcells send and receive in one direction 
    double *gcell_sendx; // left  & right
    double *gcell_sendy; // front & back
    double *gcell_sendz; // down  & up 

    double *gcell_recvx;
    double *gcell_recvy;
    double *gcell_recvz;

    
    // time slice information
    int mytid;
    int tsnum; 

    Grid(PinT* conf); 
    virtual ~Grid();

    virtual void init() {
        dummy("INIT");
    }
    // for easily initialize the grid variables
    inline void set_val4all(size_t ind, double val) {
        //u_f[ind] = val;      // for fine solver ,  not necessary 
        //u_c[ind] = val;      // for coarse solver, not necessary 
        u_cprev[ind] = val;  // for coarse solver previous time iteration  
        u_start[ind] = val;  // for start point of the current time slice
        u_end[ind] = val;    // for end point of the current time slice
    }

    // boundary condition
    void bc();
    void bc(double *d); 
    void bc_1d(double *d); 
    void bc_2d(double *d); 
    void bc_3d(double *d); 

    // create the grid topology
    void create_topology(); 
    void create_topology_1d();
    void create_topology_2d();
    void create_topology_3d();

    // guardcell exchange  
    void guardcell();
    void guardcell(double *d);
    void guardcell_1d(double *d);
    void guardcell_2d(double *d);
    void guardcell_3d(double *d);

    // grid topology. left:X:right; front:Y:back; down:Z:up
    int left, right, front, back, down, up;
  
    /**
     * get the XYZ global geographical value from their local outer index,
     * often used in iniialization and position-related calculations.
     */
    // ind is the outer index
    inline double getX(int ind, bool outer=true){
        size_t gind;
        if(outer)
            gind = (idx+ind-nguard);
        else gind = idx+ind;
        return gind*dx + dx/2; 
    }
    inline double getY(int ind, bool outer=true ){
        if(ndim<2) return 0;
        size_t gind;
        if(outer)
            gind = (idy+ind-nguard);
        else gind = idy + ind;

        return gind*dy + dy/2;
    }
    inline double getZ(int ind, bool outer=true){
        if(ndim<3) return 0;
        size_t gind; 
        if(outer)
            gind = (idz+ind-nguard);
        else gind = idz + ind;

        return gind*dz + dz/2;
    }
    inline double getGX(int igx, int ind){
        return (igx+ind+0.5)*dx; 
    }
    inline double getGY(int igy, int ind){
        if(ndim<2) return 0;
        return (igy+ind+0.5)*dy; 
    }
    inline double getGZ(int igz, int ind){
        if(ndim<3) return 0;
        return (igz+ind+0.5)*dz; 
    }

    // get the linear index of the 1-2-3D index, the set of functions are not frequently used.
    // if the ndim<3, the extra parameters will be ignored 
    // ix, iy, iz is the outer coords
    inline long getInnerIdxO(int ix, int iy, int iz){
        if(ndim==3) return ny*nx*(iz-nguard)+nx*(iy-nguard)+(ix-nguard);
        else if(ndim==2) return nx*(iy-nguard)+(ix-nguard);
        else return ix-nguard;
    }
    // ix, iy, iz is the inner coords
    inline long getInnerIdxI(int ix, int iy, int iz){
        if(ndim==3) return ny*nx*iz + nx*iy + ix;
        else if(ndim==2) return nx*iy + ix;
        else return ix;
    }
    // ix, iy, iz is the outer coords
    inline long getOuterIdx(int ix, int iy, int iz){
        if(ndim==3) return sy*sx*iz + sx*iy + ix;
        else if(ndim==2) return sx*iy + ix;
        else return ix;
    }
    
    // get rid of guard cells and pack inner grid data into a buffer  
    inline void pack_data(double *p, double *buf){
        switch(ndim) {
          case 1: pack_1d_(nxyz, &nguard, p, buf); break;
          case 2: pack_2d_(nxyz, &nguard, p, buf); break;
          case 3: pack_3d_(nxyz, &nguard, p, buf); break; 
        }
    }

    virtual void bc_1d_l(double *d){ dummy("BC"); } 
    virtual void bc_1d_r(double *d){ dummy("BC"); } 

    virtual void bc_2d_l(double *d){ dummy("BC"); } 
    virtual void bc_2d_r(double *d){ dummy("BC"); }
    virtual void bc_2d_f(double *d){ dummy("BC"); }
    virtual void bc_2d_b(double *d){ dummy("BC"); }

    virtual void bc_3d_l(double *d){ dummy("BC"); } 
    virtual void bc_3d_r(double *d){ dummy("BC"); }
    virtual void bc_3d_f(double *d){ dummy("BC"); }
    virtual void bc_3d_b(double *d){ dummy("BC"); }
    virtual void bc_3d_d(double *d){ dummy("BC"); }
    virtual void bc_3d_u(double *d){ dummy("BC"); }

    // space domain, for only one datum with double type and SUM operation 
    void sp_allreduce(double *d);  //d is input and output
    void sp_allreduce(double *d, double *o); //d is input, o is output
    // time-space domain
    void allreduce(double *d, double *o, MPI_Op op);
    
    // the basic output function, write grid local variable for debug, 
    // the wrapper of Output.var_..._Z, writing the local variables into debug file only for the last time slice  
    // only for X-Y cross sections, if for X-Z or Y-Z cross sections, turn to the Output class  
    void output_local(double *data, bool inner_only);

    void output_local_h5(double *data); // for large inner data only

    // aggregate all the final results from all the grids within the same space domain and write to disk 
    // ASCII mode only for X-Y cross sections along the Z direction   
    //
    // NOTE: 
    // For large scale performance test on HPC environment, it is best not to use the global output, 
    // the function has not yet used parallel I/O of HDF5, 
    // it will take quite some time for collecting data and writing them serially at master process.     
    void output_global(const char* flagname, bool h5=false);
    inline void output_global_h5(const char* flagname) {
        output_global(flagname, true);
    }

};
#endif
