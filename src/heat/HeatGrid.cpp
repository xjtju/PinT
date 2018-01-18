#include "HeatGrid.h"

HeatGrid::HeatGrid(PinT *conf) : Grid(conf){
} 

int HeatGrid::init() {
    if(this->ndim == 1) init1d();
    if(this->ndim == 2) init2d();
    if(this->ndim == 3) init3d();
}
void HeatGrid::init1d(){
    long ind = 0;
    double x, unk;
    for(int i = nguard; i<nx+nguard ; i++){
        ind = i;
        x = this->getX(i);  // global coordincate of the cell
        unk = cos(2*x);
        
        // set the variables used by Parareal method 
        u_f[ind] = unk;      // for fine solver  
        u_c[ind] = unk;      // for coarse solver  
        u_cprev[ind] = unk;  // for coarse solver previous time iteration  
        u_start[ind] = unk;  // for start point of the current time slice
        u_end[ind] = unk;    // for end point of the current time slice
    }
}

void HeatGrid::init2d(){
    long ind = 0;     
    double xdist, ydist, unk;
    for(int j = nguard; j<ny+nguard; j++)
    for(int i = nguard; i<nx+nguard; i++){
       // get distance from the center of the whole geographical domain
       xdist = this->getX(i) -  conf->Xspan/2 ; 
       ydist = this->getY(j) -  conf->Yspan/2 ;
       
       ind = this->getOuterIdx(i, j, 0);  // get the index including the guard cell 
        
       if( abs(xdist)<=0.2 &&  abs(ydist)<=0.2 )
           unk = 100.0; 
       else unk = 0.0; 
       
       // set the variables used by Parareal method 
       u_f[ind] = unk;      // for fine solver  
       u_c[ind] = unk;      // for coarse solver  
       u_cprev[ind] = unk;  // for coarse solver previous time iteration 
       u_start[ind] = unk;  // for start point of the current time slice
       u_end[ind] = unk;    // for end point of the current time slice
    }
}

void HeatGrid::init3d(){
    long ind = 0;     
    double xdist, ydist, zdist, unk;
    for(int k = nguard; k<nz+nguard; k++)
    for(int j = nguard; j<ny+nguard; j++)
    for(int i = nguard; i<nx+nguard; i++){
       // get distance from the center of the whole geographical domain
       xdist = this->getX(i) -  conf->Xspan/2 ; 
       ydist = this->getY(j) -  conf->Yspan/2 ;
       zdist = this->getZ(k) -  conf->Zspan/2 ;

       ind = this->getOuterIdx(i, j, k);  // get the index including the guard cell 
        
       if( abs(xdist)<=0.2 &&  abs(ydist)<=0.2 && abs(zdist)<=0.2 )
           unk = 100.0; 
       else unk = 0.0; 
        
       // set the variables used by Parareal method 
       u_f[ind] = unk;      // for fine solver  
       u_c[ind] = unk;      // for coarse solver  
       u_cprev[ind] = unk;  // for coarse solver previous time iteration 
       u_start[ind] = unk;  // for start point of the current time slice
       u_end[ind] = unk;    // for end point of the current time slice
    }
}

