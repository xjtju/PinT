
#ifndef _PinT_SOR_H_
#define _PinT_SOR_H_ 1

/**
 * the template of SOR (successive over-relaxation) method for solving Ax=B.
 *
 * NOTE:
 *  Though the calculations of SOR are much less than BiCG, the convergence rate is much slow.
 *  BiCG is recommended in practice.
 *
 *  The relaxation factor, the omega value must be carefully chosen according the real problem. 
 */

#include "common.h"
#include "LS.h"
#include "blas.h"
#include "blassor.h"

class SOR : public LS{

private:

    inline void sor2_core(double *x, double *b, double *A, int *color, double *err){
         
        if(ndim==3)      sor2_core_3d_(grid->nxyz, &nguard, x, b, A, color, &omega, err) ;
        else if(ndim==2) sor2_core_2d_(grid->nxyz, &nguard, x, b, A, color, &omega, err) ;
        else if(ndim==1) sor2_core_1d_(grid->nxyz, &nguard, x, b, A, color, &omega, err) ;
    }

public:


    bool checkCnvg = true;  //default is not check
    double omega;  //relaxation factor, alias of relaxfactor

    // the wrapper red-black successive over-relaxation (sor2_core)
    int solve(double *x, double *b, double *A);
    
    SOR(PinT* c, Grid *g):LS(c,g){ 
        omega = relaxfactor; 
    }
    
    void chk_eps(double *x, double res_nrm2, double *err); 
    
    void set_omega(double omg)   { omega     = omg; }
    void set_checkCnvg(bool ifg) { checkCnvg = ifg; }
};
#endif
