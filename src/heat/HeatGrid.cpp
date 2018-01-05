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
        int gxi = (idx+i-nguard); // global coordincate of the cell
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
    double x, y, unk;
    for(int j = nguard; j<ny+nguard; j++)
    for(int i = nguard; i<nx+nguard; i++){

        int gyi = (idy+j-nguard);
        int gxi = (idx+i-nguard);

        x = gxi*dx-nx/2*dx;
        y = gyi*dy-ny/2*dy;
        ind = j*sx + i;

       //if( x*x + y*y <= 0.09 )
       //     unk = 30.0; 
       // else unk = 10.0; 

        unk = 0;

        u_f[ind] = unk; 
        u_c[ind] = unk; 
        u_cprev[ind] = unk; 
        u_start[ind] = unk;
        u_end[ind] = unk;
    }
}
