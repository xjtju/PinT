#include "PFMGrid.h"

// du/dt = D*(du/dx)^{2} + k*u(u-1.0)*(u-0.5+beta) 
//
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
    /*
    for(int i = nguard; i<nx+nguard ; i++){
        ind = i;
        x = this->getX(i);  // global coordincate of the cell

        if(ind < this->size/2 )
            unk = 1;
        else unk = 0;
        // set the variables used by Parareal method 
        //unk = cos(2*x);
        this->set_val4all(ind, unk);
    }
    */
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
