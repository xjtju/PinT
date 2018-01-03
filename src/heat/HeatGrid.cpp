#include "HeatGrid.h"

HeatGrid::HeatGrid(PinT *conf) : Grid(conf){
} 

int HeatGrid::init() {
    if(this->ndim == 1) init1d();
    if(this->ndim == 2) init2d();
}
void HeatGrid::init1d(){
    long ind = 0;
    for(int i = nguard; i<nx+nguard ; i++){
        ind = i;
        int gxi = (ind+i-nguard); // global coordincate of the cell
        double unk = cos(2*gxi*dx);
        
        u_f[ind] = unk; 
        u_c[ind] = unk; 
        u_cprev[ind] = unk; 
        u_start[ind] = unk;
        u_end[ind] = unk;
    }
}

void HeatGrid::init2d(){
    long ind = 0;     
    for(int j = ny+nguard-1; j>=nguard ; j--)
    for(int i = nguard; i<nx+nguard ; i++){

        int gyi = (idy+j-nguard);
        int gxi = (idx+i-nguard);
        ind = j*sx + i;
        double unk = cos(2*gxi*dx)*cos(2*gyi*dy);
        
        u_f[ind] = unk; 
        u_c[ind] = unk; 
        u_cprev[ind] = unk; 
        u_start[ind] = unk;
        u_end[ind] = unk;
    }
}
