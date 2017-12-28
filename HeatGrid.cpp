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
    u_f[0] = u_f[1];
    u_f[nx+nguard] = u_f[nx];
    
    u_c[0] = u_c[1];
    u_c[nx+nguard] = u_c[nx];

    u_cprev[0] = u_cprev[1];
    u_cprev[nx+nguard] = u_cprev[nx];

    u_start[0] = u_start[1];
    u_start[nx+nguard] = u_start[nx];

    u_end[0] = u_end[1];
    u_end[nx+nguard] = u_end[nx];
}
