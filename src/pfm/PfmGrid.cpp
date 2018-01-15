#include "PFMGrid.h"

PFMGrid::PFMGrid(PinT *conf) : Grid(conf){

} 

int PFMGrid::init() {
    if(this->ndim == 1) init1d();
    if(this->ndim == 2) init2d();
    if(this->ndim == 3) init3d();
}
void PFMGrid::init1d(){
    long ind = 0;
    double x, unk;
    for(int i = nguard; i<nx+nguard ; i++){
        ind = i;
        x = this->getX(i);  // global coordincate of the cell

        unk = 0.1*sin(2*PI*x) + 0.01*cos(4*PI*x) + 0.06*sin(4*PI*x) + 0.02*cos(10*PI*x);
        // set the variables used by Parareal method 
        this->set_val4all(ind, unk);
    }
}

void PFMGrid::init2d(){
    long ind = 0;     
    double xdist, ydist, unk;
    for(int j = nguard; j<ny+nguard; j++)
    for(int i = nguard; i<nx+nguard; i++){
       // get distance from the center of the whole geographical domain
       xdist = this->getX(i) -  conf->Xspan/2 ; 
       ydist = this->getY(j) -  conf->Yspan/2 ;
       
       ind = this->getOuterIdx(i, j, 0);  // get the index including the guard cell 
        
       unk = 0.0; 
       
       this->set_val4all(ind, unk);
    }
}

void PFMGrid::init3d(){
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
        
       unk = 0.0; 
        this->set_val4all(ind, unk);
    }
}
