#include "Grid.h"

Grid::Grid(PinT *conf) {
    this->conf = conf; 
    
    this->ndim = conf->ndim;

    this->nx = conf->nx;
    this->ny = conf->ny;
    this->nz = conf->nz;

    this->nxyz[0] = conf->nx;
    this->nxyz[1] = conf->ny;
    this->nxyz[2] = conf->nz;

    this->nguard = conf->nguard;
    this->bc_type = conf->bc_type;
    this->bc_val = conf->bc_val;

    this->dx = conf->dx;
    this->dy = conf->dy;
    this->dz = conf->dz;
    
    this->sx= nx + 2*nguard;
    this->sy= ny + 2*nguard;
    this->sz= nz + 2*nguard;
    //this->sy = ny;
    //this->sz = nz;
    this->ngxyz[0] = nguard;
    this->ngxyz[1] = nguard;
    this->ngxyz[2] = nguard;

    if(ndim == 1){
        this->size = this->sx;
        this->gcsx = nguard;
        this->gcell_sendx = alloc_mem(nguard);
        this->gcell_recvx = alloc_mem(nguard);
    }
    if(ndim==2) {
        this->size = this->sx * this->sy; 
        this->gcsx = this->sy*nguard;
        this->gcsy = this->sx*nguard;
        this->gcell_sendx = alloc_mem(this->gcsx);
        this->gcell_recvx = alloc_mem(this->gcsx);
        this->gcell_sendy = alloc_mem(this->gcsy);
        this->gcell_recvy = alloc_mem(this->gcsy);
    }
    if(ndim==3) {
        this->size = sx * sy * sz;
    }

    this->sxyz[0] = sx;
    this->sxyz[1] = sy;
    this->sxyz[2] = sz;

    this->spnum = conf->spnum;
    this->spnumx = conf->spnumx;
    this->spnumy = conf->spnumy;
    this->spnumz = conf->spnumz;

    this->myid = conf->myid;
    this->mysid = conf->mysid;
    this->sp_comm = conf->sp_comm;
    
    this->mytid = conf->mytid;
    this->tsnum = conf->tsnum;

    u_f = alloc_mem(size);
    u_c = alloc_mem(size);
    u_cprev = alloc_mem(size);
    u_start = alloc_mem(size);
    u_end = alloc_mem(size);
    u=u_end;

    create_topology();

    //after the topology is determined, the coordinate can also be fixed 
    this->idx = this->coords[0]*nx;
    if(ndim>=2) this->idy = this->coords[1]*ny;
    if(ndim==3) this->idz = this->coords[2]*nz;
}

Grid:: ~Grid(){
    free_mem(u_f);
    free_mem(u_c);
    free_mem(u_cprev);
    free_mem(u_start);
    free_mem(u_end);

    if(ndim == 3) {
        free_mem(gcell_sendz);
        free_mem(gcell_recvz);
    }
    if(ndim >= 2) {
        free_mem(gcell_sendy);
        free_mem(gcell_recvy);
    }
    if(ndim >=1) {
        free_mem(gcell_sendx);
        free_mem(gcell_recvx);
    }

    u = NULL;
}

void Grid::create_topology() {
    if(ndim == 1) create_topology_1d();
    if(ndim == 2) create_topology_2d();
    if(ndim == 3) 
        printf("WARNING : not completed for 3D !");
}
/**
 * create virtual space topology
 */
void Grid::create_topology_1d(){
    int periods[1] ={0};
    int dims[1];
    int coord[1]; 

    dims[0] = spnumx; 

    MPI_Cart_create(*sp_comm, ndim, dims, periods, 1, &st_comm);
    MPI_Cart_coords(st_comm, mysid, ndim, coord);
    MPI_Cart_rank(st_comm, coord, &st_rank);

    MPI_Cart_shift(st_comm, 0, 1, &left, &right);
    this->coords = coord;
    //printf("I am %d: (%d); originally %d. topology : %d | %d | %d \n",st_rank,coord_1d[0], mysid, left, st_rank, right);
}

void Grid::create_topology_2d(){
    int periods[2] = {0,0};
    int dims[2];
    int coord[2]; 

    dims[0] = spnumx; 
    dims[1] = spnumy;

    MPI_Cart_create(*sp_comm, ndim, dims, periods, 2, &st_comm);
    MPI_Cart_coords(st_comm, mysid, ndim, coord);
    MPI_Cart_rank(st_comm, coord, &st_rank);

    MPI_Cart_shift(st_comm, 0, 1, &left, &right);
    MPI_Cart_shift(st_comm, 1, 1, &front, &back);

    this->coords = coord;
//    printf("I am %d: (%d,%d); originally %d. dim0 : %d | %d | %d \n",st_rank,coord[0],coord[1], mysid, left, st_rank, right);
//    printf("I am %d: (%d,%d); originally %d. dim1 : %d | %d | %d \n",st_rank,coord[0],coord[1], mysid, front, st_rank, back);
}

/*
 * In topology, the guard cells also include boundary cells locating the whole space domain border. 
 * so the guard cell function will call the bc function at the end of it 
 */

void Grid::guardcell(double* d) {
    if(ndim == 1) guardcell_1d(d);
    if(ndim == 2) guardcell_2d(d);
}
void Grid::bc(double* d){
    if(ndim == 1) bc_1d(d);
    if(ndim == 2) bc_2d(d);
}

void Grid::guardcell_1d(double* d) {
   MPI_Request req;
   MPI_Status stat;
   int ierr;
   
   int sg = this->gcsx;

   packgc_1d_r_(nxyz, &nguard, d, gcell_sendx);
   MPI_Sendrecv(gcell_sendx, sg, MPI_DOUBLE, right,  8008, 
                gcell_recvx, sg, MPI_DOUBLE, left, 8008, st_comm, &stat);
   unpackgc_1d_l_(nxyz, &nguard, d, gcell_recvx);

   //printf("%d L: %f, %f ", st_rank, gcell_send[0], gcell_recv[0]);
 
   packgc_1d_l_(nxyz, &nguard, d, gcell_sendx);
   MPI_Sendrecv(gcell_sendx, sg, MPI_DOUBLE, left, 9009, 
                gcell_recvx, sg, MPI_DOUBLE, right,  9009, st_comm, &stat);
   unpackgc_1d_r_(nxyz, &nguard, d, gcell_recvx);

   bc_1d(d);
}

void Grid::guardcell_2d(double* d) {
   MPI_Request req;
   MPI_Status stat;
   int sg;
 
   // X: left - right : sy * nguard  
   sg = this->gcsx; 

   if(MPI_PROC_NULL!=right)  // not right border
       packgc_2d_r_(nxyz, &nguard, d, gcell_sendx);
   MPI_Sendrecv(gcell_sendx, sg, MPI_DOUBLE, right,  8008, 
                gcell_recvx, sg, MPI_DOUBLE, left, 8008, st_comm, &stat);
   if(MPI_PROC_NULL!=left)  // not left border
       unpackgc_2d_l_(nxyz,&nguard,d,gcell_recvx);
    
   if(MPI_PROC_NULL!=left)  // not left border
       packgc_2d_l_(nxyz, &nguard, d, gcell_sendx);
   MPI_Sendrecv(gcell_sendx, sg, MPI_DOUBLE, left, 9009, 
                gcell_recvx, sg, MPI_DOUBLE, right,  9009, st_comm, &stat);
   if(MPI_PROC_NULL!=right)  // not right border
       unpackgc_2d_r_(nxyz, &nguard, d, gcell_recvx);

   // Y: front - back  
   sg = this->gcsy;
   
   if(MPI_PROC_NULL!=back)  
       packgc_2d_b_(nxyz, &nguard, d, gcell_sendy);
   MPI_Sendrecv(gcell_sendy, sg, MPI_DOUBLE, back,  6006, 
                gcell_recvy, sg, MPI_DOUBLE, front, 6006, st_comm, &stat);
   if(MPI_PROC_NULL!=front)  
       unpackgc_2d_f_(nxyz,&nguard,d,gcell_recvy);
    
   if(MPI_PROC_NULL!=front)  
       packgc_2d_f_(nxyz, &nguard, d, gcell_sendy);
   MPI_Sendrecv(gcell_sendy, sg, MPI_DOUBLE, front, 7007, 
                gcell_recvy, sg, MPI_DOUBLE, back,  7007, st_comm, &stat);
   if(MPI_PROC_NULL!=back)  
       unpackgc_2d_b_(nxyz, &nguard, d, gcell_recvy);

   //printf("sg=%d, size=%d, ng=%d, nx=%d, ny=%d, gsbx=%f.\n", sg,size,nguard,nxyz[0], nxyz[1], gcell_sendx[sg-1]);
   bc_2d(d);
}

/**
 * Deprecated.
 * Usually it is not necessary to synchonize guard cells for all the variables.
 */
void Grid::guardcell() {
    guardcell(u_f);
    guardcell(u_c);
    guardcell(u_cprev);
    guardcell(u_start);
    guardcell(u_end);
}
/**
 * Deprecated.
 * Usually it is not necessary to set border cells for all the variables,
 * although the bc function is very light-weight (no MPI communication).
 */
void Grid::bc(){
    bc(u_f); 
    bc(u_c);   
    bc(u_cprev);   
    bc(u_start);   
    bc(u_end);    
}
// nguard = 1, now space parallel is not considered
void Grid::bc_1d(double* d) {
    int ng = nguard;
    if( 0==bc_type ){ //fixed value
       if(MPI_PROC_NULL==left)  // left border
           bc_val_1d_l_(nxyz, &ng, d, &bc_val);
       if(MPI_PROC_NULL==right)  //right border
           bc_val_1d_r_(nxyz, &ng, d, &bc_val);
    }
    else if( 1==bc_type ) { //reflected
       if(MPI_PROC_NULL==left)  
           bc_ref_1d_l_(nxyz, &ng, d);
       if(MPI_PROC_NULL==right)  
           bc_ref_1d_r_(nxyz, &ng, d);
   }
}
void Grid::bc_2d(double* d) {
    int ng = nguard;
    if( 0==bc_type ) {
       if(MPI_PROC_NULL==left)   bc_val_2d_l_(nxyz, &ng, d, &bc_val);
       if(MPI_PROC_NULL==right)  bc_val_2d_r_(nxyz, &ng, d, &bc_val);
       if(MPI_PROC_NULL==front)  bc_val_2d_f_(nxyz, &ng, d, &bc_val);
       if(MPI_PROC_NULL==back)   bc_val_2d_b_(nxyz, &ng, d, &bc_val);
    }
    else if( 1==bc_type ) { //reflected
       if(MPI_PROC_NULL==left)   bc_ref_2d_l_(nxyz, &ng, d);
       if(MPI_PROC_NULL==right)  bc_ref_2d_r_(nxyz, &ng, d);
       if(MPI_PROC_NULL==front)  bc_ref_2d_f_(nxyz, &ng, d);
       if(MPI_PROC_NULL==back)   bc_ref_2d_b_(nxyz, &ng, d);
    }
}
/**
 * MPI_Allreduce double
 */
void Grid::sp_allreduce(double *d) {
    double tmp[1]; 
    MPI_Allreduce(d, tmp, 1, MPI_DOUBLE, MPI_SUM, *sp_comm);
    *d = *tmp;
}

void Grid::sp_allreduce(double *d, double *o) {
    MPI_Allreduce(d, o, 1, MPI_DOUBLE, MPI_SUM,*sp_comm);
}

void Grid::allreduce(double *d, double *o, int op) {
    MPI_Allreduce(d, o, 1, MPI_DOUBLE, op, MPI_COMM_WORLD); 
}

/**
 * collect the final result from all space domain and output 
 * it simply write the u_end to file, used only for debug ,not for massive running.
 * file name : mytid.mysid.txt 
 */
void Grid::output() {
    if (mytid != (tsnum-1)) return;  //only output the last time slice

    FILE * fp;
    char fname[20];
    int ind = 0;   
    sprintf(fname, "%s_%d.%d.txt", conf->debug_pre,mytid,mysid); 

    fp = fopen (fname,"w");
    for(int j=sy-1; j>=0 ; j--) { 
        if( (j==nguard-1) || (j==sy-nguard-1)) fprintf(fp, "  ------------  \n");
        for(int i = 0; i < sx ; i++){
            ind = j*sx + i;
            if( (i==nguard) || (i==sx-nguard)) fprintf(fp, " | ");
            fprintf (fp, "  %12.8f  ", u_end[ind]);
        }
        fprintf(fp,"\n"); 
    }
    fprintf(fp,"\n"); 

    fclose (fp);
}
