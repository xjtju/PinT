
#ifndef _PinT_SOR_H_
#define _PinT_SOR_H_ 1

/**
 * the template of SOR (successive over-relaxation) method for solving Ax=B.
 *
 * NOTE:
 *  Though the calculations of SOR are much less than BiCG, the convergence rate is much slow.
 *  BiCG is recommended in practice.
 *
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
    // these control parammeters can be over-writen by sub classes
    double eps = 1.0e-6;
    // SOR : relaxation factor, the omega value must be carefully chosen according the real problem. 
    // for heat equation, the optimal value may be 1.7,
    // but for Allen-Cahn equation, 1.0 may be a better choice.
    double  omega= 1.0;      

    int itmax = 20; 
    bool checkCnvg = false;  //default is not check

    // the wrapper red-black successive over-relaxation (sor2_core)
    int solve(double *x, double *b, double *A);
    
    SOR(PinT* c, Grid *g):LS(c,g){ }
    
    void chk_eps(double *x, double res_nrm2, double *err); 

    void set_omega(double omg) {
        omega = omg;
    }

    void set_eps(double err) {
        eps = err;
    }
    void set_itmax(int iter) {
        itmax = iter;
    }
    void set_checkCnvg(bool ifg) {
        checkCnvg = ifg;
    }
};
#endif
