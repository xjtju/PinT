#include "Grid.h"

Grid::Grid(int nx, int nguard, double dx, double dt) {
    this->nx = nx;
    this->nguard = nguard;
    this->dx = dx;
    this->dt = dt;
    this->size = nx+2*nguard;
}

Grid:: ~Grid(){
    free_mem(x);
}
