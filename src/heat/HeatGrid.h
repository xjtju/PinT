#ifndef PinT_HEATGRID_H_
#define PinT_HEATGRID_H_ 1 

#include "Grid.h"

/**
 * the 1D Heat Diffusion with reflect boundary condition
 *   Ut = kUxx, k = .061644
 */
class HeatGrid : public Grid {

public: 

    HeatGrid(PinT *conf); 
    int init();
};
#endif
