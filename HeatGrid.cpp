#include "HeatGrid.h"

HeatGrid::HeatGrid(int nx, int ng, double dx, double dt) : Grid(nx,ng,dx,dt){} 

int HeatGrid::init(){
    x = alloc_mem(size);
         
    for(int i = nguard; i<nx+nguard ; i++){
       x[i] = cos(2*i*dx); 
    }
    bc();
    return 0;
}

// nguard = 1
void HeatGrid::bc(){
    x[0] = x[1];
    x[nx+nguard] = x[nx];
}
