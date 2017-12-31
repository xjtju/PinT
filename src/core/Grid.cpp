#include "Grid.h"

Grid::Grid(PinT *conf) {
    this->conf = conf; 
    
    this->dims = conf->dims;

    this->nx = conf->nx;
    this->ny = conf->ny;
    this->nz = conf->nz;
    this->nxyz[0] = conf->nx;
    this->nxyz[1] = conf->ny;
    this->nxyz[2] = conf->nz;

    this->nguard = conf->nguard;

    this->dx = conf->dx;
    this->dy = conf->dy;
    this->dz = conf->dz;
    
    this->sx= nx + 2*nguard;
    this->sy = ny;
    this->sz = nz;
    this->ngxyz[0] = nguard;
    this->ngxyz[1] = 0;
    this->ngxyz[2] = 0;
    if(dims>=2) {
        this->ngxyz[1] = nguard;
        this->sy= ny + 2*nguard;
    }
    if(dims==3) {
        this->ngxyz[2] = nguard;
        this->sz= nz + 2*nguard;
    }

    
    this->sxyz[0] = sx;
    this->sxyz[1] = sy;
    this->sxyz[2] = sz;
    this->size = sx * sy * sz;

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
