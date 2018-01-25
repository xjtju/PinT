#include "BD4Solver.h"


void BD4Solver::init_holder() {
    soln_2 = alloc_mem(this->size);
    soln_3 = alloc_mem(this->size);
    soln_4 = alloc_mem(this->size);
}
