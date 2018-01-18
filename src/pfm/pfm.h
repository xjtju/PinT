#ifndef _PinT_PFM_H_
#define _PinT_PFM_H_ 1

extern "C" {

    void rhs_ac_1d_(int *nxyz, double *lamda, int *ng, double *b, double *soln, double *soln_, double *g1, 
            double *theta, double *dtk, double *beta_);

    void rhs_g1_ac_1d_(int *nxyz, double *lamda, int *ng, double *soln, double *g1, 
            double *theta, double *dtk, double *beta_);

    void stencil_ac_1d_(int *nxyz, double *lamda, int *ng, double *bcp, double *soln,
            double *theta, double *dtk, double *beta_);

    void update_ac_1d_(int *nxyz, int *ng, double *soln, double *delta);

    void bc_ac_1d_(int *nxyz, int *ng, double *soln);

}
#endif
