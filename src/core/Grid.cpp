#include "Grid.h"
#include "Output.h"
/**
 * X (0 ... i) : left          right 
 * Y (0 ... j) : front         back 
 * Z (0 ... k) : down(under) up     
 */
Grid::Grid(PinT *conf) {
    this->conf = conf; 
    
    this->ndim = conf->ndim;

    this->nguard = conf->nguard;
    this->bc_type = conf->bc_type;
    this->bc_val = conf->bc_val;

    this->nx = conf->nx;
    this->ny = conf->ny;
    this->nz = conf->nz;

    this->nxyz[0] = conf->nx;
    this->dx = conf->dx;
    this->sx= nx + 2*nguard;
    this->ngxyz[0] = nguard;
    this->sy = 1; // TEST
    this->sz = 1; // TEST
    if(ndim>=2) { 
        this->nxyz[1] = conf->ny;
        this->dy = conf->dy;
        this->sy= ny + 2*nguard;
        this->ngxyz[1] = nguard;
    }
    if(ndim==3) {
        this->nxyz[2] = conf->nz;
        this->dz = conf->dz;
        this->sz= nz + 2*nguard;
        this->ngxyz[2] = nguard;
    }

    if(ndim == 1){
        this->inner_size = nx;
        this->size = this->sx;
        this->gcsx = nguard;
        this->gcell_sendx = alloc_mem(nguard);
        this->gcell_recvx = alloc_mem(nguard);
    }
    if(ndim==2) {
        this->inner_size = nx*ny;
        this->size = this->sx * this->sy; 
        this->gcsx = this->sy*nguard;
        this->gcsy = this->sx*nguard;
        this->gcell_sendx = alloc_mem(this->gcsx);
        this->gcell_recvx = alloc_mem(this->gcsx);
        this->gcell_sendy = alloc_mem(this->gcsy);
        this->gcell_recvy = alloc_mem(this->gcsy);
    }
    if(ndim==3) {
        this->inner_size = nx*ny*nz;
        this->size = sx * sy * sz;
        this->gcsx = this->sz*this->sy*nguard;
        this->gcsy = this->sz*this->sx*nguard;
        this->gcsz = this->sy*this->sx*nguard;
        this->gcell_sendx = alloc_mem(this->gcsx);
        this->gcell_recvx = alloc_mem(this->gcsx);
        this->gcell_sendy = alloc_mem(this->gcsy);
        this->gcell_recvy = alloc_mem(this->gcsy);
        this->gcell_sendz = alloc_mem(this->gcsz);
        this->gcell_recvz = alloc_mem(this->gcsz);
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

    this->coords = new int[ndim]; 
    this->coords_= new int[ndim]; 
    this->periods= new int[ndim];
    this->dims   = new int[ndim];  

    create_topology();

    // after the topology is determined, the coordinate can also be fixed 
    this->idx = this->coords[0]*nx;
    if(ndim>=2) this->idy = this->coords[1]*ny;
    if(ndim==3) this->idz = this->coords[2]*nz;
    
    // allocate memory for grid level variables 
    u_f = alloc_mem(size);
    u_c = alloc_mem(size);
    u_cprev = alloc_mem(size);
    u_start = alloc_mem(size);
    u_end = alloc_mem(size);
    u=u_end;
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
    if(ndim == 3) create_topology_3d();
}
/**
 * create virtual space topology
 */
void Grid::create_topology_1d(){
    periods[0] = 0;
    dims[0] = spnumx; 

    MPI_Cart_create(*sp_comm, ndim, dims, periods, 1, &st_comm);
    MPI_Cart_coords(st_comm, mysid, ndim, coords);
    MPI_Cart_rank(st_comm, coords, &st_rank);

    MPI_Cart_shift(st_comm, 0, 1, &left, &right);
}

void Grid::create_topology_2d(){
    periods[0] = 0;
    periods[1] = 0;
    dims[0] = spnumx; 
    dims[1] = spnumy;

    MPI_Cart_create(*sp_comm, ndim, dims, periods, 2, &st_comm);
    MPI_Cart_coords(st_comm, mysid, ndim, coords);
    MPI_Cart_rank(st_comm, coords, &st_rank);

    MPI_Cart_shift(st_comm, 0, 1, &left,  &right);
    MPI_Cart_shift(st_comm, 1, 1, &front, &back);
}

void Grid::create_topology_3d(){
    periods[0] = 0;
    periods[1] = 0;
    periods[2] = 0;
    dims[0] = spnumx; 
    dims[1] = spnumy;
    dims[2] = spnumz;

    MPI_Cart_create(*sp_comm, ndim, dims, periods, 3, &st_comm);
    MPI_Cart_coords(st_comm, mysid, ndim, coords);
    MPI_Cart_rank(st_comm, coords, &st_rank);

    MPI_Cart_shift(st_comm, 0, 1, &left,   &right);
    MPI_Cart_shift(st_comm, 1, 1, &front,  &back);
    MPI_Cart_shift(st_comm, 2, 1, &down, &up);

    //printf("I am %d: (%d, %d, %d); originally %d. X : %d | %d | %d \n",st_rank,coords[0],coords[1], coords[2], mysid, left, st_rank, right);
    //printf("I am %d: (%d, %d, %d); originally %d. Y : %d | %d | %d \n",st_rank,coords[0],coords[1], coords[2], mysid, front, st_rank, back);
    //printf("I am %d: (%d, %d, %d); originally %d. Z : %d | %d | %d \n",st_rank,coords[0],coords[1], coords[2], mysid, down, st_rank, up);
}
/*
 * In topology, the guard cells also include boundary cells locating the whole space domain border. 
 * But the guard cell function will not automatically call the bc function at the end of it,
 * because boundary condition is usually applied between two time steps, but guard cell may be exchanged within one time steps when necessary. 
 */
void Grid::guardcell(double* d) {
    //monitor.start(Monitor::GC);

    if(ndim == 1) guardcell_1d(d);
    else if(ndim == 2) guardcell_2d(d);
    else if(ndim == 3) guardcell_3d(d);

    //monitor.stop(Monitor::GC);
}

void Grid::guardcell_1d(double* d) {
   MPI_Request req;
   MPI_Status stat;
   int ierr;
   
   int sg = this->gcsx;

   if(MPI_PROC_NULL!=right)  // not right border
       packgc_1d_r_(nxyz, &nguard, d, gcell_sendx);
   MPI_Sendrecv(gcell_sendx, sg, MPI_DOUBLE, right,  8008, 
                gcell_recvx, sg, MPI_DOUBLE, left, 8008, st_comm, &stat);
   if(MPI_PROC_NULL!=left)  // not left border
      unpackgc_1d_l_(nxyz, &nguard, d, gcell_recvx);

   //printf("%d L: %f, %f ", st_rank, gcell_send[0], gcell_recv[0]);
 
   if(MPI_PROC_NULL!=left)  // not left border
       packgc_1d_l_(nxyz, &nguard, d, gcell_sendx);
   MPI_Sendrecv(gcell_sendx, sg, MPI_DOUBLE, left, 9009, 
                gcell_recvx, sg, MPI_DOUBLE, right,  9009, st_comm, &stat);
   if(MPI_PROC_NULL!=right)  // not right border
       unpackgc_1d_r_(nxyz, &nguard, d, gcell_recvx);
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
}

//
// WARN : 3D guardcell uses sxyz not nxyz !!
//
void Grid::guardcell_3d(double* d) {
   MPI_Request req;
   MPI_Status stat;
   int sg;
 
   // X: left - right : Y-Z surface  
   sg = this->gcsx; 
   if(MPI_PROC_NULL!=right)  // not right border
       packgc_3d_r_(sxyz, &nguard, d, gcell_sendx);
   MPI_Sendrecv(gcell_sendx, sg, MPI_DOUBLE, right,  8008, 
                gcell_recvx, sg, MPI_DOUBLE, left, 8008, st_comm, &stat);
   if(MPI_PROC_NULL!=left)  // not left border
       unpackgc_3d_l_(sxyz,&nguard,d,gcell_recvx);
    
   if(MPI_PROC_NULL!=left)  // not left border
       packgc_3d_l_(sxyz, &nguard, d, gcell_sendx);
   MPI_Sendrecv(gcell_sendx, sg, MPI_DOUBLE, left, 9009, 
                gcell_recvx, sg, MPI_DOUBLE, right,  9009, st_comm, &stat);
   if(MPI_PROC_NULL!=right)  // not right border
       unpackgc_3d_r_(sxyz, &nguard, d, gcell_recvx);

   
   // Y: front - back : X-Z surface 
   sg = this->gcsy;
   if(MPI_PROC_NULL!=back)  
       packgc_3d_b_(sxyz, &nguard, d, gcell_sendy);
   MPI_Sendrecv(gcell_sendy, sg, MPI_DOUBLE, back,  6006, 
                gcell_recvy, sg, MPI_DOUBLE, front, 6006, st_comm, &stat);
   if(MPI_PROC_NULL!=front)  
       unpackgc_3d_f_(sxyz,&nguard,d,gcell_recvy);
    
   if(MPI_PROC_NULL!=front)  
       packgc_3d_f_(sxyz, &nguard, d, gcell_sendy);
   MPI_Sendrecv(gcell_sendy, sg, MPI_DOUBLE, front, 7007, 
                gcell_recvy, sg, MPI_DOUBLE, back,  7007, st_comm, &stat);
   if(MPI_PROC_NULL!=back)  
       unpackgc_3d_b_(sxyz, &nguard, d, gcell_recvy);

    
   // Z: down - up : X-Y surface 
   sg = this->gcsz;
   if(MPI_PROC_NULL!=up)  
       packgc_3d_u_(sxyz, &nguard, d, gcell_sendz);
   MPI_Sendrecv(gcell_sendz, sg, MPI_DOUBLE, up,  4004, 
                gcell_recvz, sg, MPI_DOUBLE, down, 4004, st_comm, &stat);
   if(MPI_PROC_NULL!=down)  
       unpackgc_3d_d_(sxyz,&nguard,d,gcell_recvz);
    
   if(MPI_PROC_NULL!=down)  
       packgc_3d_d_(sxyz, &nguard, d, gcell_sendz);
   MPI_Sendrecv(gcell_sendz, sg, MPI_DOUBLE, down, 5005, 
                gcell_recvz, sg, MPI_DOUBLE, up,  5005, st_comm, &stat);
   if(MPI_PROC_NULL!=up)  
       unpackgc_3d_u_(sxyz, &nguard, d, gcell_recvz);
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
void Grid::bc(double* d){
    if(ndim == 1) bc_1d(d);
    if(ndim == 2) bc_2d(d);
    if(ndim == 3) bc_3d(d);
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

void Grid::bc_3d(double *d){
    int ng = nguard;
    if( 0==bc_type) {
        if(MPI_PROC_NULL==left)   bc_val_3d_l_(sxyz, &ng, d, &bc_val);
        if(MPI_PROC_NULL==right)  bc_val_3d_r_(sxyz, &ng, d, &bc_val);
        if(MPI_PROC_NULL==front)  bc_val_3d_f_(sxyz, &ng, d, &bc_val);
        if(MPI_PROC_NULL==back)   bc_val_3d_b_(sxyz, &ng, d, &bc_val);
        if(MPI_PROC_NULL==down)   bc_val_3d_d_(sxyz, &ng, d, &bc_val);
        if(MPI_PROC_NULL==up)     bc_val_3d_u_(sxyz, &ng, d, &bc_val);
    }
    else if( 1==bc_type ){
        if(MPI_PROC_NULL==left)   bc_ref_3d_l_(sxyz, &ng, d);
        if(MPI_PROC_NULL==right)  bc_ref_3d_r_(sxyz, &ng, d);
        if(MPI_PROC_NULL==front)  bc_ref_3d_f_(sxyz, &ng, d);
        if(MPI_PROC_NULL==back)   bc_ref_3d_b_(sxyz, &ng, d);
        if(MPI_PROC_NULL==down)   bc_ref_3d_d_(sxyz, &ng, d);
        if(MPI_PROC_NULL==up)     bc_ref_3d_u_(sxyz, &ng, d);
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
void Grid::output_local(double *p, bool inner_only) {
    if (mytid != (tsnum-1)) return;  //only output the last time slice

    Output out = Output(this);

    FILE * fp;
    char fname[30];
    int ind = 0;   
    double *buf;

    memset(fname, 0, sizeof(char)*30);
    sprintf(fname, "%s_%d.%d.txt", conf->debug_pre,mytid,mysid); 

    fp = fopen (fname,"w");
    
    if(inner_only){
        buf = alloc_mem(this->inner_size); 
        pack_data(p, buf);  // get rid of the guard cell
        out.var_inner_Z(fp, buf, true);
        free_mem(buf);
    }else { 
        out.var_outer_Z(fp, p);
    }

    fclose (fp);
}

/**
 * at current, the function is not well tuned, but it can work.
 */
void Grid::output_local_h5(double *p) {
    if (mytid != (tsnum-1)) return;  //only output the last time slice

    Output out = Output(this);

    FILE * fp;
    char fname[30];
    int ind = 0;   
    double *buf;

    memset(fname, 0, sizeof(char)*30);
    sprintf(fname, "%s_%d.%d.h5", conf->debug_pre,mytid,mysid); 
        
    buf = alloc_mem(this->inner_size); 
    pack_data(p, buf); 
    out.write_h5(fname, buf); 
    free_mem(buf);
}

// file name : mytid.all.txt 
void Grid::output_global(){
    if (mytid != (tsnum-1)) return;  //only output the last time slice
    
    MPI_Request req, req_p;
    MPI_Status sta, sta_p;
    int ierr ;
    int source, dest;

    FILE *fp;
    char fname[21];
    int ind = 0;   
    
    int idbuf[3];
    double* sendrecv_buf = alloc_mem(this->inner_size); 

    dest = spnum - 1;
    if( mysid != dest ){ //is not the last space grid 
        printf("[%d.%d] is sending data to the last grid [%d.%d]\n", mytid, mysid,mytid, dest);
        pack_data(u_end, sendrecv_buf);
        ierr = MPI_Isend(sendrecv_buf, inner_size, MPI_DOUBLE, dest, 9999, *sp_comm, &req); 
        MPI_Wait(&req, &sta);
        return;
    }

    // the output task is left to the last space grid 
    sprintf(fname, "%s_%d.all.txt", conf->debug_pre,mytid); 
    fp = fopen (fname,"w");
    Output out = Output(this);
    for(int k=0; k<spnumz; k++)
    for(int j=0; j<spnumy; j++)
    for(int i=0; i<spnumx; i++) {
        source = i + j*spnumx + k*spnumy*spnumx;
        if(source == mysid) continue;
        printf("[%d.%d] is aggregating data from the front grid [%d.%d]\n", mytid, mysid, mytid, source);
        MPI_Irecv(sendrecv_buf, inner_size, MPI_DOUBLE, source, 9999, *sp_comm, &req_p);
        MPI_Wait(&req_p, &sta_p);

        MPI_Cart_coords(st_comm, source, ndim, coords_);
        out.coord(fp, coords_);
        out.var_inner_Z(fp,sendrecv_buf, true);
    }
    // output itself 
    pack_data(u_end, sendrecv_buf); // get rid of guard cell
    out.coord(fp, coords);
    out.var_inner_Z(fp,sendrecv_buf,true);

    fclose(fp);

    free_mem(sendrecv_buf);

    printf("All the data has been aggregated into the file : %s \n", fname);
}

