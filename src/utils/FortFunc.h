#ifndef _FORTRAN_FUNC_H_
#define _FORTRAN_FUNC_H_ 1

extern "C" {
    //
    // 1D : 2 directions
    //
    void bc_val_1d_l_(int* nxyz, int *ng, double* p, double *val);
    void bc_val_1d_r_(int* nxyz, int *ng, double* p, double *val);

    void bc_ref_1d_l_(int* nxyz, int *ng, double* p);
    void bc_ref_1d_r_(int* nxyz, int *ng, double* p);

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

    void bc_ref_2d_l_(int* nxyz, int *ng, double* p);
    void bc_ref_2d_r_(int* nxyz, int *ng, double* p);
    void bc_ref_2d_f_(int* nxyz, int *ng, double* p);
    void bc_ref_2d_b_(int* nxyz, int *ng, double* p);

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
    // 3D : 6 directions
    //

    // left  right 
    // front back
    // top   bottom/under
    //
    
    // for pack and unpack all inner grid data for aggregating result 

    void pack_1d_(int* nxyz, int *ng, double* p, double* gdata);

    void pack_2d_(int* nxyz, int *ng, double* p, double* gdata);
}
#endif

