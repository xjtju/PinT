#include "Grid.h"

Grid::Grid(int nx, int ng, double dx, double dt) {
    this->nx = nx;
    this->nguard = ng;
    this->dx = dx;
    this->dt = dt;
    this->size = nx+2*nguard;
}

Grid:: ~Grid(){
    free_mem(x);
}
