#ifndef _PinT_PFM_H_
#define _PinT_PFM_H_ 1

extern "C" {
    //
    // For Crank-Nicolson and Newton-Raphson
    //
    // 1D
    void rhs_ac_1d_(int *nxyz, double *lamdaxyz, int *ng, double *b, double *soln, double *soln_, double *g1, 
            double *theta, double *dtk, double *beta_);

    void rhs_g1_ac_1d_(int *nxyz, double *lamdaxyz, int *ng, double *soln, double *g1, 
            double *theta, double *dtk, double *beta_);

    void stencil_ac_1d_(int *nxyz, double *lamdaxyz, int *ng, double *A, double *soln,
            double *theta, double *dtk, double *beta_);


    // customized border condition
    void bc_pfm_ac_1d_l_(int* nxyz, int *ng, double* soln);  // left border 
    void bc_pfm_ac_1d_r_(int* nxyz, int *ng, double* soln);  // right border

    // 2D
    void rhs_ac_2d_(int *nxyz, double *lamdaxyz, int *ng, double *b, double *soln, double *soln_, double *g1, 
            double *theta, double *dtk, double *beta_);

    void rhs_g1_ac_2d_(int *nxyz, double *lamdaxyz, int *ng, double *soln, double *g1, 
            double *theta, double *dtk, double *beta_);

    void stencil_ac_2d_(int *nxyz, double *lamdaxyz, int *ng, double *A, double *soln,
            double *theta, double *dtk, double *beta_);


    // only for example, not used 
    void bc_pfm_ac_2d_l_(int* sxyz, int *ng, double* soln);  // left  border : west 
    void bc_pfm_ac_2d_r_(int* sxyz, int *ng, double* soln);  // right border : east 
    void bc_pfm_ac_2d_f_(int* sxyz, int *ng, double* soln);  // front border : south 
    void bc_pfm_ac_2d_b_(int* sxyz, int *ng, double* soln);  // back  border : north

    //
    // 3D
    //
    void rhs_ac_3d_(int *nxyz, double *lamdaxyz, int *ng, double *b, double *soln, double *soln_, double *g1, 
            double *theta, double *dtk, double *beta_);

    void rhs_g1_ac_3d_(int *nxyz, double *lamdaxyz, int *ng, double *soln, double *g1, 
            double *theta, double *dtk, double *beta_);

    void stencil_ac_3d_(int *nxyz, double *lamdaxyz, int *ng, double *A, double *soln,
            double *theta, double *dtk, double *beta_);


    //
    //For first order forward Euler 
    //
    void euler_rhs_ac_3d_(int *nxyz, double *lamdaxyz, int *ng, double *b, double *soln, double *dtk, double *beta_);
    void euler_rhs_ac_2d_(int *nxyz, double *lamdaxyz, int *ng, double *b, double *soln, double *dtk, double *beta_);
    void euler_rhs_ac_1d_(int *nxyz, double *lamdaxyz, int *ng, double *b, double *soln, double *dtk, double *beta_);
    

    //
    // For 4th order backward differentiation method of AC
    //
    // 1D
    void rhs_g1_ac_bd4_1d_(int *nxyz, int *ng, double *soln_1, double *soln_2, double *soln_3, double *soln_4, 
            double *g1, double *dtk, double *beta_);
    
    void rhs_ac_bd4_1d_(int *nxyz, double *lamdaxyz, int *ng, double *b, double *soln, double *soln_, 
            double *g1, double *dtk, double *beta_);

    void stencil_ac_bd4_1d_(int *nxyz, double *lamdaxyz, int *ng, double *A, double *soln, double *dtk, double *beta_);

    // 2D
    void rhs_g1_ac_bd4_2d_(int *nxyz, int *ng, double *soln_1, double *soln_2, double *soln_3, double *soln_4, 
            double *g1, double *dtk, double *beta_);
    
    void rhs_ac_bd4_2d_(int *nxyz, double *lamdaxyz, int *ng, double *b, double *soln, double *soln_, 
            double *g1, double *dtk, double *beta_);

    void stencil_ac_bd4_2d_(int *nxyz, double *lamdaxyz, int *ng, double *A, double *soln, double *dtk, double *beta_);

    // 3D
    void rhs_g1_ac_bd4_3d_(int *nxyz, int *ng, double *soln_1, double *soln_2, double *soln_3, double *soln_4, 
            double *g1, double *dtk, double *beta_);
    
    void rhs_ac_bd4_3d_(int *nxyz, double *lamdaxyz, int *ng, double *b, double *soln, double *soln_, 
            double *g1, double *dtk, double *beta_);

    void stencil_ac_bd4_3d_(int *nxyz, double *lamdaxyz, int *ng, double *A, double *soln, double *dtk, double *beta_);


}
#endif
