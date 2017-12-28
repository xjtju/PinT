#include "Grid.h"

Grid::Grid(PinT *conf) {
    this->conf = conf; 
    this->nx = conf->sub_nx;
    this->nguard = conf->nguard;
    this->dx = conf->dx;
    this->size = nx+2*nguard;

    u_f = alloc_mem(size);
    u_c = alloc_mem(size);
    u_cprev = alloc_mem(size);
    u_start = alloc_mem(size);
    u_end = alloc_mem(size);
    u=u_end;
}

Grid:: ~Grid(){
    free_mem(u_f);
    free_mem(u_c);
    free_mem(u_cprev);
    free_mem(u_start);
    free_mem(u_end);
    u = NULL;
}
