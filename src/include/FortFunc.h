#ifndef _FORTRAN_FUNC_H_
#define _FORTRAN_FUNC_H_ 1

extern "C" {
    //
    // 1D : 2 directions
    //
    void bc_val_1d_l_(int* nxyz, int *ng, double* p, double *val);
    void bc_val_1d_r_(int* nxyz, int *ng, double* p, double *val);

    void bc_der_1d_l_(int* nxyz, int *ng, double* p);
    void bc_der_1d_r_(int* nxyz, int *ng, double* p);

    void packgc_1d_l_(int* nxyz, int *ng, double* p, double* gdata);
    void packgc_1d_r_(int* nxyz, int *ng, double* p, double* gdata);

    void unpackgc_1d_l_(int* nxyz, int *ng, double* p, double* gdata);
    void unpackgc_1d_r_(int* nxyz, int *ng, double* p, double* gdata);

    // 
    // 2D : 4 directions
    //  
    void bc_val_2d_l_(int* nxyz, int *ng, double* p, double *val);
    void bc_val_2d_r_(int* nxyz, int *ng, double* p, double *val);
    void bc_val_2d_f_(int* nxyz, int *ng, double* p, double *val);
    void bc_val_2d_b_(int* nxyz, int *ng, double* p, double *val);

    void bc_der_2d_l_(int* nxyz, int *ng, double* p);
    void bc_der_2d_r_(int* nxyz, int *ng, double* p);
    void bc_der_2d_f_(int* nxyz, int *ng, double* p);
    void bc_der_2d_b_(int* nxyz, int *ng, double* p);

    // left right : X
    void packgc_2d_l_(int* nxyz, int *ng, double* p, double* gdata);
    void packgc_2d_r_(int* nxyz, int *ng, double* p, double* gdata);

    void unpackgc_2d_l_(int* nxyz, int *ng, double* p, double* gdata);
    void unpackgc_2d_r_(int* nxyz, int *ng, double* p, double* gdata);

    // front back : Y
    void packgc_2d_f_(int* nxyz, int *ng, double* p, double* gdata);
    void packgc_2d_b_(int* nxyz, int *ng, double* p, double* gdata);

    void unpackgc_2d_f_(int* nxyz, int *ng, double* p, double* gdata);
    void unpackgc_2d_b_(int* nxyz, int *ng, double* p, double* gdata);


    // 
    // 3D : 6 directions. left<->right; front<->back; down<->up 
    // NOTE : sxyz (outer size) is used, not the nxyz (inner size) 
    //
    // left right : X
    void packgc_3d_l_(int* sxyz, int *ng, double* p, double* gdata);
    void packgc_3d_r_(int* sxyz, int *ng, double* p, double* gdata);
    void unpackgc_3d_l_(int* sxyz, int *ng, double* p, double* gdata);
    void unpackgc_3d_r_(int* sxyz, int *ng, double* p, double* gdata);

    void bc_val_3d_l_(int* sxyz, int *ng, double* p, double *val);
    void bc_val_3d_r_(int* sxyz, int *ng, double* p, double *val);
    void bc_der_3d_l_(int* sxyz, int *ng, double* p );
    void bc_der_3d_r_(int* sxyz, int *ng, double* p );

    // front back : Y
    void packgc_3d_f_(int* sxyz, int *ng, double* p, double* gdata);
    void packgc_3d_b_(int* sxyz, int *ng, double* p, double* gdata);
    void unpackgc_3d_f_(int* sxyz, int *ng, double* p, double* gdata);
    void unpackgc_3d_b_(int* sxyz, int *ng, double* p, double* gdata);

    void bc_val_3d_f_(int* sxyz, int *ng, double* p, double *val);
    void bc_val_3d_b_(int* sxyz, int *ng, double* p, double *val);
    void bc_der_3d_f_(int* sxyz, int *ng, double* p );
    void bc_der_3d_b_(int* sxyz, int *ng, double* p );
  

    // down up : Z
    void packgc_3d_d_(int* sxyz, int *ng, double* p, double* gdata);
    void packgc_3d_u_(int* sxyz, int *ng, double* p, double* gdata);
    void unpackgc_3d_d_(int* sxyz, int *ng, double* p, double* gdata);
    void unpackgc_3d_u_(int* sxyz, int *ng, double* p, double* gdata);
   

    void bc_val_3d_d_(int* sxyz, int *ng, double* p, double *val);
    void bc_val_3d_u_(int* sxyz, int *ng, double* p, double *val);
    void bc_der_3d_d_(int* sxyz, int *ng, double* p );
    void bc_der_3d_u_(int* sxyz, int *ng, double* p );



    // for pack and unpack all inner grid data for aggregating result 
    void pack_1d_(int* nxyz, int *ng, double* p, double* gdata);
    void pack_2d_(int* nxyz, int *ng, double* p, double* gdata);
    void pack_3d_(int* nxyz, int *ng, double* p, double* gdata);

    // common update function for soln[n+1] = soln[n] + delta
    void update_soln_1d_(int *nxyz, int *ng, double *soln, double *delta);
    void update_soln_3d_(int *nxyz, int *ng, double *soln, double *delta);
    void update_soln_2d_(int *nxyz, int *ng, double *soln, double *delta);

}
#endif

