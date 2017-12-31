#include "HeatGrid.h"

HeatGrid::HeatGrid(PinT *conf) : Grid(conf){
} 

int HeatGrid::init(){
         
    for(int i = nguard; i<nx+nguard ; i++){
        double unk = cos(2*i*dx);
        u_f[i] = unk; 
        u_c[i] = unk; 
        u_cprev[i] = unk; 
        u_start[i] = unk;
        u_end[i] = unk;
    }
    bc();
    return 0;
}

// nguard = 1, now space parallel is not considered
void HeatGrid::bc(){
    bc_(nxyz, ngxyz, u_f); 
    bc_(nxyz, ngxyz, u_c);   
    bc_(nxyz, ngxyz, u_cprev);   
    bc_(nxyz, ngxyz, u_start);   
    bc_(nxyz, ngxyz, u_end);    
}
