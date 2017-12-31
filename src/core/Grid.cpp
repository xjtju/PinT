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
    this->mysid = conf->myid;
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
    printf("gcell size : %d\n", sguard);
}

Grid:: ~Grid(){
    free_mem(u_f);
    free_mem(u_c);
    free_mem(u_cprev);
    free_mem(u_start);
    free_mem(u_end);
    u = NULL;
}
void Grid::create_topology(){
    int dims[1];
    int periods[1];
    int coord_1d[1]; 

    dims[0] = spnumx; 
    periods[0] = 0;

    MPI_Cart_create(*sp_comm, ndim, dims, periods, 1, &comm1d);
    MPI_Cart_coords(comm1d, mysid, ndim, coord_1d);
    MPI_Cart_rank(comm1d, coord_1d, &rank_1d);

    printf("I am %d: (%d); originally %d\n",rank_1d,coord_1d[0],mysid);
}
void Grid::guardcell() {
   MPI_Request req;
   MPI_Status stat;
   int ierr;
   int left, right, front, back, top, bottom;

   MPI_Cart_shift(comm1d,0,-1,&rank_1d,&left);
   MPI_Cart_shift(comm1d,0,+1,&rank_1d,&right);

   MPI_Sendrecv(gcell_send, sguard, MPI_DOUBLE, left,  8008, 
                gcell_recv, sguard, MPI_DOUBLE, right, 8008, comm1d, &stat);
   MPI_Sendrecv(gcell_send, sguard, MPI_DOUBLE, right, 9009, 
                gcell_recv, sguard, MPI_DOUBLE, left,  9009, comm1d, &stat);

}
