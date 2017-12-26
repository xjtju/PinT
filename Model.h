
#ifndef PinT_MODEL_H
#define PinT_MODEL_H 1

#include "common.h"

/**
 * in current, the model is only for holding the physical variables.
 * In the further, it can be extended to a mesh structure
 */

struct Model {

public:
    int nx;
    int nguard;
    int size;

    double *x;

    Model(int nx, int nguard) {
        this->nx = nx;
        this->nguard = nguard;
        this->size = nx+2*nguard;

        x = alloc_mem(size);
    }

    ~Model(){
        free_mem(x);
    }
     
    virtual int init_x()=0;

    //boundary condition
    virtual void bc()=0;
};

#endif
