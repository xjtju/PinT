#ifndef PinT_PFMGRID_H_
#define PinT_PFMGRID_H_ 1 

#include "Grid.h"
#include "pfm.h"
#include "PFMParams.h"

/**
 *  
 *  Phase Field Model Grid
 *
 *  Initialize solution variables for 1D/2D/3D
 * 
 *  The class is also an example of customized BC
 *  specific boundary condition, when the default bc functions of Grid cannot satisfy the requirement. 
 */
class PFMGrid : public Grid {
   
    PFMParams param;

public: 

    PFMGrid(PinT *conf):Grid(conf){
    }
    virtual ~PFMGrid() {};

    void init();
    void init1d();
    void init2d();
    void init3d();

    // specific boundary condition
    void bc_1d_l(double *d){ bc_pfm_ac_1d_l_(nxyz, &nguard, d); }
    void bc_1d_r(double *d){ bc_pfm_ac_1d_r_(nxyz, &nguard, d); }
   
    // NOTE: 2d is not tested 
    void bc_2d_l(double *d){ bc_pfm_ac_2d_l_(sxyz, &nguard, d); }
    void bc_2d_r(double *d){ bc_pfm_ac_2d_l_(sxyz, &nguard, d); }
    void bc_2d_f(double *d){ bc_pfm_ac_2d_f_(sxyz, &nguard, d); }
    void bc_2d_b(double *d){ bc_pfm_ac_2d_b_(sxyz, &nguard, d); }

};
#endif
