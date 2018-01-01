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
    this->sy = ny;
    this->sz = nz;
    this->ngxyz[0] = nguard;
    this->ngxyz[1] = 0;
    this->ngxyz[2] = 0;
    if(ndim>=2) {
        this->ngxyz[1] = nguard;
        this->sy= ny + 2*nguard;
    }
    if(ndim==3) {
        this->ngxyz[2] = nguard;
        this->sz= nz + 2*nguard;
    }

    this->sxyz[0] = sx;
    this->sxyz[1] = sy;
    this->sxyz[2] = sz;
    this->size = sx * sy * sz;

    this->spnum = conf->spnum;
    this->spnumx = conf->spnumx;
    this->spnumy = conf->spnumy;
    this->spnumz = conf->spnumz;

    this->myid = conf->myid;
    this->mysid = conf->mysid;
    this->sp_comm = conf->sp_comm;

    u_f = alloc_mem(size);
    u_c = alloc_mem(size);
    u_cprev = alloc_mem(size);
    u_start = alloc_mem(size);
    u_end = alloc_mem(size);
    u=u_end;

    sguard = nguard * sy * sz;  
    gcell_send = alloc_mem(sguard);
    gcell_recv = alloc_mem(sguard);

    //printf("gcell size : %d\n", sguard);
    create_topology();

    //after the topology is determined, the coordinate can also be fixed 
    this->idx = rank_1d*nx;
}

Grid:: ~Grid(){
    free_mem(u_f);
    free_mem(u_c);
    free_mem(u_cprev);
    free_mem(u_start);
    free_mem(u_end);
    u = NULL;
}

/**
 * create virtual space topology
 */
void Grid::create_topology(){
    int dims[1];
    int periods[1];
    int coord_1d[1]; 

    dims[0] = spnumx; 
    periods[0] = 0;

    MPI_Cart_create(*sp_comm, ndim, dims, periods, 1, &comm1d);
    MPI_Cart_coords(comm1d, mysid, ndim, coord_1d);
    MPI_Cart_rank(comm1d, coord_1d, &rank_1d);

    MPI_Cart_shift(comm1d, 0, 1, &left, &right);
    printf("I am %d: (%d); originally %d. topology : %d | %d | %d \n",rank_1d,coord_1d[0], mysid, left, rank_1d, right);
}

void Grid::guardcell(double* d) {
   MPI_Request req;
   MPI_Status stat;
   int ierr;

   gcell_send[0] = d[sx-2];
   MPI_Sendrecv(gcell_send, sguard, MPI_DOUBLE, right,  8008, 
                gcell_recv, sguard, MPI_DOUBLE, left, 8008, comm1d, &stat);
   d[0] = gcell_recv[0];

   //printf("%d L: %f, %f ", rank_1d, gcell_send[0], gcell_recv[0]);
 
   gcell_send[0] = d[1];
   MPI_Sendrecv(gcell_send, sguard, MPI_DOUBLE, left, 9009, 
                gcell_recv, sguard, MPI_DOUBLE, right,  9009, comm1d, &stat);
   d[sx-1] = gcell_recv[0];

   bc(d);
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

// nguard = 1, now space parallel is not considered
void Grid::bc(double* d){
   if( 0==bc_type ){
       //fixed value
   }else ( 1==bc_type ) {
       //reflected
   }

   if(MPI_PROC_NULL==left)  // I am is the most left border
      d[0] = d[1]; 
   if(MPI_PROC_NULL==right)  //right border
      d[sx-1] = d[sx-2]; 
}

void Grid::bc(){
    bc(u_f); 
    bc(u_c);   
    bc(u_cprev);   
    bc(u_start);   
    bc(u_end);    
}
/**
 * collect the final result from all space domain and output 
 */
void Grid::output() {
    
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
