#include "HeatGrid.h"

HeatGrid::HeatGrid(PinT *conf) : Grid(conf){
} 

int HeatGrid::init(){
         
    for(int i = nguard; i<nx+nguard ; i++){
        int gxi = (idx+i-nguard); // global coordincate of the cell
        double unk = cos(2*gxi*dx);
        u_f[i] = unk; 
        u_c[i] = unk; 
        u_cprev[i] = unk; 
        u_start[i] = unk;
        u_end[i] = unk;
    }
    //guardcell();
    return 0;
}

