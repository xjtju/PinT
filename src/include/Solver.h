#ifndef PinT_SOLVER_H_
#define PinT_SOLVER_H_ 

#include "PinT.h"
#include "Grid.h"
#include "LS.h" 
#include "SOR.h"
#include "PBiCGStab.h"

/**
 * the abstract interface of fine/coarse solvers for PinT framework
 * the interface is designed as simple as possible to adapt to more problem-specific calcaluations.
 *
 * mainly designed for COARSE solver
 *
 * For correctly transfering the solutions between adjacent time slices for (non)standard parareal, 
 * Some particular virtual functions are designed,
 * the basic class only provides default implementations for standard parareal algorithm.
 * 
 * For non-standard parareal example, see pfm/BD4 solver  
 *
 */
class Solver {
private:
    void init(PinT *conf, Grid *g) {
        this->conf = conf; 
        this->grid = g;

        this->ndim = g->ndim;

        this->nx = g->nx;
        this->ny = g->ny;
        this->nz = g->nz;

        this->sx = g->sx;
        this->sy = g->sy;
        this->sz = g->sz;

        this->nguard = g->nguard;

        this->inner_size = g->inner_size; 
        this->outer_size = g->size;

        b  = alloc_mem(this->size);
        if(ndim==1) A = alloc_mem(3*this->size);  // 3-point stencil for 1D
        if(ndim==2) A = alloc_mem(5*this->size); 
        if(ndim==3) A = alloc_mem(7*this->size); 
     }

protected:
    Grid *grid;
    PinT *conf;
    
    double *b;   // RHS, b of Ax=b
    double *A;   // stencil matrix A of Ax=b 

    // all the following variables are the same with the corresponding one in the Grid
    // holding the most being used variables here is just for convenience only 
    int ndim;

    int nx;
    int ny;
    int nz;

    int sx;
    int sy;
    int sz;

    int nguard;
    
    size_t inner_size; // guard cells not included 
    size_t outer_size; // guard cells included
    size_t &size = outer_size;      // alias of outer_size 

    int steps;     // the number of time steps in one time slice

    bool isFine = true;  // is fine solver or coarse solver, default is true

    LS *hypre = NULL;  // linear solver 
   
    /**
     * sub classes can overwrite it, create any proper linear solver
     *
     * SOR  has less calculations but slow convergence rate;
     * BiCG has more calculations but fast convergence rate.
     *
     * NOTE: 
     *   It is not necessary for sub classes to free the LS instance. 
     */
    virtual LS* getLS(PinT *conf, Grid *grid) {
        if(PinT::LS_SOR_ID  == conf->ls_solver){
            if(grid->myid==0) printf("INFO: default linear solver is SOR \n");
            return new SOR(conf, grid);
        }
        else if(PinT::LS_BiCG_ID == conf->ls_solver ) 
        {
            if(grid->myid==0) printf("INFO: default linear solver is BiCG \n");
            return new PBiCGStab(conf, grid); 
        }
        else { 
            if(grid->myid==0) fprintf(stderr, "WARN : default linear solver is not set \n");
            return NULL; 
        }
    } 

public:
    // set a proper linear solver 
    inline void set_LS(LS *ls){
        hypre = ls;
    }

    Solver(PinT *conf, Grid *g){
        init(conf, g);
    }

    Solver(PinT *conf, Grid *g, bool isFS){
        init(conf, g);

        // set the steps for time integrating (evolving)
        this->isFine = isFS;
        if(isFine)
            this->steps = conf->f_steps;  
        else this->steps = conf->c_steps;
    }
    virtual ~Solver() {
        free_mem(b);
        free_mem(A);

        delete hypre; 

        if(grid->myid==0 && conf->verbose)
        printf("INFO: The memory allocated by the base solver has been released.\n");
    }
     
     // get the current solution 
     virtual double* getSoln() {
         if(isFine) return grid->u_f;
         else return grid->u_c;
     }
     
     // set the initial value
     virtual void init() {
         printf("WARN: the blank init function is used for solver");
     } 

     // integrate over one time slice, that is one time slice iteration (each Kpar)
     // return total iteration count
     virtual unsigned long  evolve() = 0;

     /***************** non-standard parareal extension interfaces **********************/

     // return the solution of previous iteration of coarse solver 
     virtual double* prev_solns() {
         if(isFine) return NULL;
         return grid->u_cprev;
     }

     virtual double* curr_solns(){
         return getSoln(); 
     }

     // preserve the solution of coarse solver
     virtual void backup_prevs(){
        if(isFine) return ;     
        blas_cp_(grid->u_cprev, grid->u_c, &size);  
     }

     virtual size_t solnsize() { // the size of send/recv variables
         return size;
     }
    
     // return the pointer of variables to be sended    
     virtual double* sendslns() {
         return grid->u_end; 
     }
     // return the pointer of variables to be received 
     virtual double* recvslns() {
         return grid->u_start;  
     }
     // pack and unpack the send/recv data, nothing need to be done in the standard parareal 
     // It is not an easy task to make the pair functions general, not completed
     virtual void pack()   { }
     virtual void unpack() { }
     // write back the final solution (u_end) from correction item (Driver.pint_sum) 
     virtual void update_uend(){ }
};

#endif
