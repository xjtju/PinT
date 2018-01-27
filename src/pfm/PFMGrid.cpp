#include "PFMGrid.h"

// set the initial value
void PFMGrid::init() {

    param.init(conf,this);
    if(myid==0)
        param.print();

    if(ndim == 1) init1d();
    if(ndim == 2) init2d();
    if(ndim == 3) init3d();
}

// set the initial value
void PFMGrid::init1d() {
    double d = param.d;
    double k = param.k;

    double xi = sqrt(2.0*d/k);
    double val;
    double midx = conf->Xspan/2 ; 
    double xdist;
    for(int i=nguard; i<nx+nguard; i++){
        double x = this->getX(i); 
        val = 1.0 + tanh( (x-0.5)/xi);    // initial value 
        val = 1.0 - 0.5*val;
        //xdist = x - midx; 
        //if(abs(xdist) < 0.2) 
        //   val = 0.8;
        //else val = 0;
        this->set_val4all(i,val);
    }
}

void PFMGrid::init2d(){
    long ind = 0;     
    double xdist, ydist, unk;
    for(int j = nguard; j<ny+nguard; j++)
    for(int i = nguard; i<nx+nguard; i++){
       xdist = this->getX(i) -  conf->Xspan/2 ; 
       ydist = this->getY(j) -  conf->Yspan/2 ;
       
       ind = this->getOuterIdx(i, j, 0);  
       if( abs(xdist)<=0.4 &&  abs(ydist)<=0.4 )
           unk = 1.0; 
       else unk = 0.0; 
       
       this->set_val4all(ind, unk);
    }
}

void PFMGrid::init3d(){
    long ind = 0;     
    double xdist, ydist, zdist, unk;
    for(int k = nguard; k<nz+nguard; k++)
    for(int j = nguard; j<ny+nguard; j++)
    for(int i = nguard; i<nx+nguard; i++){
       xdist = this->getX(i) -  conf->Xspan/2 ; 
       ydist = this->getY(j) -  conf->Yspan/2 ;
       zdist = this->getZ(k) -  conf->Zspan/2 ;
       ind = this->getOuterIdx(i, j, k);  
       if( (abs(xdist)<=0.4) &&  (abs(ydist)<=0.4) && (abs(zdist)<=0.4) )
           unk = 1.0; 
       else unk = 0.0; 
       
       this->set_val4all(ind, unk);
    }
}
/**
 * problem specific parameter parsing
 */
int pfm_inih(void* obj, const char* section, const char* name, const char* value) {

    PFMParams* pfm = (PFMParams*)obj;
    if (MATCH("pfm","ac_kval"))        { pfm->k = atof(value); }
    if (MATCH("pfm","ac_dval"))        { pfm->d = atof(value); }
    if (MATCH("pfm","ac_beta"))        { pfm->beta = atof(value); }
    if (MATCH("pfm", "newton_itmax"))  { pfm->newton_itmax = atoi(value); }
    if (MATCH("pfm", "theta"))         { pfm->theta = atof(value); }
    if (MATCH("pfm", "csolver"))         { pfm->csolver = atof(value); }
    if (MATCH("pfm", "fsolver"))         { pfm->fsolver = atof(value); }

    return 0;
}

