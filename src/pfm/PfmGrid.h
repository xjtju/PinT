#ifndef PinT_PFMGRID_H_
#define PinT_PFMGRID_H_ 1 

#include "Grid.h"

/**
 * Phase Field Model, Allen-Cahnn equation  
 */
class PFMGrid : public Grid {

public: 

    PFMGrid(PinT *conf); 
    int init();
    void init1d();
    void init2d();
    void init3d();
};
#endif
