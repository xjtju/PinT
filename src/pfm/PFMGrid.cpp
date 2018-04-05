#include "PFMGrid.h"

/**
 * NOTE:  
 *  1. The final (steady) shape of result is related with the initial value.
 *
 *  2. The steady condition is different for 1D/2D/3D due to the diffuse effect and reactive energy item. 
 *     So the final results are also different at 1D/2D/3D even if the program is absolutely correct.
 *     See PFMParams.h for details.
 */

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

    double r = 0.4;
    double xi = sqrt(2.0*d/k);
    double val;
    double midx = conf->Xspan/2 ; 
    double xdist;
    for(int i=nguard; i<nx+nguard; i++){
        double x = this->getX(i); 
//        val = 1.0 + tanh( (x-0.5)/xi);    // initial value 
//        val = 1.0 - 0.5*val;
        xdist = this->getX(i) -  conf->Xspan/2 ; 
        if(fabs(xdist) < r) 
           val = 1.0;
        else val = 0;
        this->set_val4all(i,val);
    }
}

void PFMGrid::init2d(){
    long ind = 0;     
    double r = 0.4;
    double xdist, ydist, unk;
    for(int j = nguard; j<ny+nguard; j++)
    for(int i = nguard; i<nx+nguard; i++){
       xdist = this->getX(i) -  conf->Xspan/2 ; 
       ydist = this->getY(j) -  conf->Yspan/2 ;
       
       ind = this->getOuterIdx(i, j, 0);  
       if( fabs(xdist)<=r ) // && fabs(ydist)<=r )
       //if( xdist<=0 )
           unk = 1.0; 
       else unk = 0.0; 
       
       this->set_val4all(ind, unk);
    }
}

void PFMGrid::init3d(){
    long ind = 0;     
    double xdist, ydist, zdist, unk;
    double r = 0.4;

    for(int k = nguard; k<nz+nguard; k++)
    for(int j = nguard; j<ny+nguard; j++)
    for(int i = nguard; i<nx+nguard; i++){
       xdist = this->getX(i) -  conf->Xspan/2 ; 
       ydist = this->getY(j) -  conf->Yspan/2 ;
       zdist = this->getZ(k) -  conf->Zspan/2 ;
       ind = this->getOuterIdx(i, j, k);  
       //if( (fabs(xdist)<=r)) //  && (fabs(ydist)<=r) && (fabs(zdist)<=r) )
       if(xdist<=0.0)
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

