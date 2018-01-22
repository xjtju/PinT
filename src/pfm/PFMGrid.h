#ifndef PinT_PFMGRID_H_
#define PinT_PFMGRID_H_ 1 

#include "Grid.h"
#include "pfm.h"
/**
 *  The example of customized BC
 *  specific boundary condition, the default bc functions of Grid cannot satisfy the requirement. 
 */
class PFMGrid : public Grid {

public: 

    PFMGrid(PinT *conf):Grid(conf){
    }
    virtual ~PFMGrid() {};

    void bc_1d_l(double *d){ bc_pfm_ac_1d_l_(nxyz, &nguard, d); }
    void bc_1d_r(double *d){ bc_pfm_ac_1d_r_(nxyz, &nguard, d); }
   

    void bc_2d_l(double *d){ bc_pfm_ac_2d_l_(sxyz, &nguard, d); }
    void bc_2d_l(double *d){ bc_pfm_ac_2d_l_(sxyz, &nguard, d); }
    void bc_2d_f(double *d){ bc_pfm_ac_2d_f_(sxyz, &nguard, d); }
    void bc_2d_b(double *d){ bc_pfm_ac_2d_b_(sxyz, &nguard, d); }

};
#endif
