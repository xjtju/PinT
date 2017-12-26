#ifndef PinT_HEATGRID_H
#define PinT_HEATGRID_H 1 

#include "Grid.h"

/**
 * the 1D Heat Diffusion with reflect boundary condition
 *   Ut = kUxx, k = .061644
 */
class HeatGrid : public Grid {

public: 

    HeatGrid(int nx, int ng, double dx, double dt); 
    int init();
    void bc();
};
#endif
