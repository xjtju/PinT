#ifndef _FORTRAN_FUNC_H_
#define _FORTRAN_FUNC_H_ 1

extern "C" {
    // 1D
    void bc_val_l_(int* nxyz, int *ng, double* p, double *val);
    void bc_val_r_(int* nxyz, int *ng, double* p, double *val);

    void bc_ref_l_(int* nxyz, int *ng, double* p);
    void bc_ref_r_(int* nxyz, int *ng, double* p);

    void packgc_1d_l_(int* nxyz, int *ng, double* p, double* gdata);
    void packgc_1d_r_(int* nxyz, int *ng, double* p, double* gdata);

    void unpackgc_1d_l_(int* nxyz, int *ng, double* p, double* gdata);
    void unpackgc_1d_r_(int* nxyz, int *ng, double* p, double* gdata);

    // 2D
    void bc_ref_2d_l_(int* nxyz, int *ng, double* p);
    void bc_ref_2d_r_(int* nxyz, int *ng, double* p);

    void packgc_2d_l_(int* nxyz, int *ng, double* p, double* gdata);
    void packgc_2d_r_(int* nxyz, int *ng, double* p, double* gdata);

    void unpackgc_2d_l_(int* nxyz, int *ng, double* p, double* gdata);
    void unpackgc_2d_r_(int* nxyz, int *ng, double* p, double* gdata);
}
#endif

