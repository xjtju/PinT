#include "Output.h"
#include "common.h"

void Output::var_outer(FILE *fp, double *p) {

}

// in the current version, the global idz is not supported, 
// in order to avoid confusing, when outputing the global data, it is better to set withIdz=false
// if set the withIdz=true, make sure to call output.coord() first to set the global z index correctly
void Output::var_inner_Z(FILE* fp, double *p, bool withIdz) {
    int i,j,k, ind, k_;
    if(withIdz) k_ = idz_;
    else k_ = idz;
    for(int k=0; k < nz; k++) { 
        if(with_coord) fprintf(fp, " z global index : %d \n", k_ + k );
        for(int j=ny-1; j>=0 ; j--) { 
            if(with_coord) fprintf(fp, "%4d :", idy_+j);  // print global id of Y
            for(int i = 0; i < nx ; i++){
                ind = grid->getInnerIdxI(i, j, k); 
                fprintf (fp," %10.5f  ", p[ind]);
            }
            fprintf(fp, "\n");
        }
        if(with_coord) for(int i = 0; i < nx ; i++) fprintf(fp, " %11d ", idx_+i);  //print global id of X
        fprintf(fp,"\n\n\n");
    }
}

/**
 * 3D
 * it is hard to express 3D data to an ASCII file, so
 * the func outputs all cross surfaces(section) orderly from one direction.   
 * direction: 0-X (the YZ surface), 1-Y, 2-Z
 */
void Output::var_outer_X(FILE *fp, double *p){
    int ind;
    for(int i=sx-1; i>=0; i--) {
        if( (i==nguard-1) || (i==sx-nguard-1)) fprintf(fp, "  ===========================  \n");
        fprintf(fp, " x global index : %d \n", idx + i - nguard);
        for(int k=sz-1; k>=0 ; k--) { 
            if( (k==nguard-1) || (k==sz-nguard-1)) fprintf(fp, "  ----------  \n");
            fprintf(fp, "%4d :", idz+k);   // global id of Z
            for(int j = 0; j < sy ; j++){
                ind = grid->getOuterIdx(i, j, k); 
                if( (j==nguard) || (j==sy-nguard)) fprintf(fp, " | ");
                fprintf (fp, "  %10.5f  ", p[ind]);
            }
            fprintf(fp,"\n"); 
        }
        for(int j = 0; j < sy ; j++) fprintf(fp, " %13d ", idy+j);  // global id of Y 
        fprintf(fp,"\n"); 
    }
}

void Output::var_outer_Y(FILE *fp, double *p){
    int ind;
    for(int j=sy-1; j>=0; j--) {
        if( (j==nguard-1) || (j==sy-nguard-1)) fprintf(fp, "  ===========================  \n");
        fprintf(fp, " y global index : %d \n", idy + j - nguard);
        for(int k=sz-1; k>=0 ; k--) { 
            if( (k==nguard-1) || (k==sz-nguard-1)) fprintf(fp, "  ----------  \n");
            fprintf(fp, "%4d :", idz+k);   // global id of Z 
            for(int i = 0; i < sx ; i++){
                ind = grid->getOuterIdx(i, j, k); 
                if( (i==nguard) || (i==sx-nguard)) fprintf(fp, " | ");
                fprintf (fp, "  %10.5f  ", p[ind]);
            }
            fprintf(fp,"\n"); 
        }
        for(int i = 0; i < sx ; i++) fprintf(fp, " %13d ", idx+i);  // global id of X
        fprintf(fp,"\n\n"); 
    }
}

void Output::var_outer_Z(FILE *fp, double *p){
    int ind;
    for(int k=sz-1; k>=0; k--) {
        if( (k==nguard-1) || (k==sz-nguard-1)) fprintf(fp, "  ===========================  \n");
        fprintf(fp, " z global index : %d \n", idz + k - nguard);
        for(int j=sy-1; j>=0 ; j--) { 
            if( (j==nguard-1) || (j==sy-nguard-1)) fprintf(fp, "  ----------  \n");
            fprintf(fp, "%4d :", idy+j);  // print global id of Y
            for(int i = 0; i < sx ; i++){
                ind = grid->getOuterIdx(i, j, k); 
                if( (i==nguard) || (i==sx-nguard)) fprintf(fp, " | ");
                fprintf (fp, "  %10.5f  ", p[ind]);
            }
            fprintf(fp,"\n");
        }
        for(int i = 0; i < sx ; i++) fprintf(fp, " %13d ", idx+i);  //print global id of X
        fprintf(fp,"\n\n"); 
    }
}

int Output::write_h5(char *fname, double *p) {
    int status = 0;

#ifdef _HDF5_
    hid_t file, dataset, dataspace, datatype;

    file = H5Fcreate( fname, H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT );
    datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
     
    // write solution data
    hsize_t dimsf[1];
    dimsf[0] = grid->inner_size;
    dataspace = H5Screate_simple(1, dimsf, NULL);
    dataset = H5Dcreate(file, "solution", datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, p);

    status = H5Sclose(dataspace);
    status = H5Dclose(dataset);
     
    // write coords
    hsize_t dimsc[2];
    dimsc[0] = (hsize_t)grid->inner_size; 
    dimsc[1] = (hsize_t)3;
    double *coords = alloc_mem(dimsc[0]*dimsc[1]); 
    size_t ind ;
    for(int k=0; k<nz; k++)
        for(int j=0; j<ny; j++)
            for(int i=0; i<nx; i++) {
                ind = grid->getInnerIdxI(i, j, k);
                coords[ind*3]   = grid->getX(i, false);
                coords[ind*3+1] = grid->getY(j, false);
                coords[ind*3+2] = grid->getZ(k, false);
            }
    dataspace = H5Screate_simple(2, dimsc, NULL);
    dataset = H5Dcreate(file, "coords", datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coords);
     
    status = H5Sclose(dataspace);
    status = H5Dclose(dataset);
    
    status = H5Fclose(file);

    free_mem(coords);
#else
    printf("INFO : HDF5 function is not yet activated !\n");
#endif   
    return status; 
}
