#include "SOR.h"
#include "Driver.h"
/**
 * the wrapper red-black successive over-relaxation (sor2_core)
 *
 * check convergency flag is not used. if used, the previous solution must be hold for comparing. 
 */
int SOR::solve(double *x, double *b, double *A) {
    bool ifg = false; 
    double res_nrm2 = 0.0;
    double err = 0.0;
    int lc = 1;
    for(lc=1; lc<=itmax; lc++) {
        for(int color=0; color<2; color++) {
            if(lc>1) grid->guardcell(x);
            sor2_core(x, b, A, &color, &res_nrm2);
         }

        if(checkCnvg) {
            chk_eps(x, res_nrm2, &err);
            if( err < eps) { 
                ifg = true;
                break;
            }
        }
    }
    if(checkCnvg && (!ifg) && force_abort) Driver::Abort("SOR loop does not converge, eps=%e\n" , err);
    return lc;
}

void SOR::chk_eps(double *x, double res_nrm2, double *err) {
    double x_nrm2 = 1.0;
    if(ndim==3) { 
        blas_dot_3d_(grid->nxyz, &nguard, x, x, &x_nrm2 );
    } else if(ndim==2) {
        blas_dot_2d_(grid->nxyz, &nguard, x, x, &x_nrm2 );
    } else if(ndim==1) {
        blas_dot_1d_(grid->nxyz, &nguard, x, x, &x_nrm2 );
    }

    grid->sp_allreduce(&res_nrm2);
    grid->sp_allreduce(&x_nrm2);

    *err = sqrt(res_nrm2/x_nrm2);
}
