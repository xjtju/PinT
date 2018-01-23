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

int Output::open_h5(char *fname){
#ifdef _HDF5_
    gfile = H5Fcreate( fname, H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT );
    if(gfile < 0) { 
        fprintf(stderr, "ERROR : create global h5 file failure\n");
        return -1;
    }
    
    gdatatype = H5Tcopy(H5T_NATIVE_DOUBLE);
    hsize_t dimsoln[1], dimcoord[2];

    dimsoln[0] = (hsize_t)grid->inner_size*spnum;  // all grids 
    hsize_t dataspace = H5Screate_simple(1, dimsoln, NULL);
    gdataset_d = H5Dcreate(gfile, "solution", gdatatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    gdimsoln[0]  = (hsize_t) grid->inner_size;  // single grid
    gmemspace_d = H5Screate_simple(1, gdimsoln,  NULL);   

   

    dimcoord[0] = (hsize_t)grid->inner_size*spnum; // all grids
    dimcoord[1] = (hsize_t)3;
    dataspace   = H5Screate_simple(2, dimcoord, NULL);
    gdataset_c = H5Dcreate(gfile, "coords", gdatatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    gdimcoord[0] = (hsize_t) grid->inner_size;  // single grid
    gdimcoord[1] = (hsize_t)3;
    gmemspace_c = H5Screate_simple(2, gdimcoord, NULL);   



    coords_ = alloc_mem(gdimcoord[0]*gdimcoord[1]);
    gridcounter = 0;
   
    //printf("solution size: %d, coords size: %d \n", dimcoord[0], dimcoord[0]*dimcoord[1]);
#else
    printf("INFO : HDF5 function is not yet activated !\n");
#endif
    return 0;
}

/**
 * NOTE : it must be used with to3dcoord(int *) pairly for calculating the coordinates correctly.
 */
int Output::append_h5(double *p){
    int status = 0;
#ifdef _HDF5_
    hid_t dataspace_id;
    
    if(gridcounter>=spnum) {
        fprintf(stderr, "ERROR : data is more than the pre-allocated size !\n");
        return -1;
    }
    //append solution variables 
    hsize_t offset[1] = {0};
    offset[0] = grid->inner_size*gridcounter;

    dataspace_id = H5Dget_space(gdataset_d);
    status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, gdimsoln, NULL);
    status = H5Dwrite (gdataset_d, H5T_NATIVE_DOUBLE, gmemspace_d, dataspace_id, H5P_DEFAULT, p);    

    H5Sclose(dataspace_id);
    
    //append coordinates 
    size_t ind ;
    for(int k=0; k<nz; k++)
        for(int j=0; j<ny; j++)
            for(int i=0; i<nx; i++) {
                ind = grid->getInnerIdxI(i, j, k);
                coords_[ind*3]   = grid->getGX(idx_, i); 
                coords_[ind*3+1] = grid->getGY(idy_, j);
                coords_[ind*3+2] = grid->getGZ(idz_, k);
            }
    hsize_t offset_[2] = {0, 0};
    offset_[0] = offset[0];
    offset_[1] = 0; 
    //printf("offset (%d , %d) \n", offset_[0], offset_[1]);
    dataspace_id = H5Dget_space(gdataset_c);
    status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset_, NULL, gdimcoord, NULL);
    status = H5Dwrite (gdataset_c, H5T_NATIVE_DOUBLE, gmemspace_c, dataspace_id, H5P_DEFAULT, coords_);    
    
    H5Sclose(dataspace_id);
    gridcounter ++; // for the next appending
#endif
    return status;
}

int Output::close_h5() {
    int status = 0;
#ifdef _HDF5_
    
    free_mem(coords_);

    status = H5Sclose(gmemspace_d);
    status = H5Sclose(gmemspace_c);
    status = H5Dclose(gdataset_c);
    status = H5Dclose(gdataset_d);
    
    status = H5Fclose(gfile); 
#endif
    return status;
}
