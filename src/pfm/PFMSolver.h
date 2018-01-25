#ifndef _PinT_PFMSOLVER_H_
#define _PinT_PFMSOLVER_H_ 1

#include "PBiCGStab.h"
#include "pfm.h"

/**
 * Phase Field Model using Newton-Raphson and Crank-Nicolson
 *
 * Allen-Cahn Equation (AC)
 *
 * Unlike the classic Heat equation with constant diffuse coefficient, 
 * the discrete formula of the AC equation based on Crank-Nicolson is nonlinear, 
 * so Newton-Raphson method is used to linearize the numerical formula of AC, and the linear solver can be applied to AC.
 * Compared to the HEAT, it is introduced extra three variables for data structure preservation during the calculation 
 * due to its nonlinear feature.
 *
 * If using the default configuration, k=16000, d=1, beta=-0.128, 
 * only after 0.1 second, the system has already reached its steady state.
 * 
 */

class PFMSolver : public Solver {

protected:
    void setup();

    virtual LS* getLS(PinT *conf, Grid *grid){
        return new PBiCGStab(conf, grid);
    }

public:

    double k;     // interfacial width related
    double d;     // diffuse coefficient 
    double beta;  // potential engrgy parameter
    double beta_;       // 0.5-beta

    double theta = 0.5;  //Crank-Nicolson; 0: Ex.Euler; 1: Im.Euler
    int newton_itmax = 5;

    double dtk;         // dt*k
    double lamda_x;     // d*dt/(dx**2)
    double lamda_y = 0.0;     
    double lamda_z = 0.0;    
    double lamdaxyz[3]; 

    double *soln;   // the current solution, pointer to the grid->u_f/u_c
    // structure preservation 
    double *soln_;  // the holder of -F^{k-1} in Newton's method when applying to nonlinear systems of equations 
    double *G1;     // the partial RHS of Newton's method for PFM  
    double *unk;    // the unknown 'x' of Ax=b, that is (Xn+1 - Xn) for Newton's method 

    //double ls_eps;
    //double ls_itmax;

    LS *hypre;  // linear solver 

    PFMSolver(PinT *c, Grid *g); 
    PFMSolver(PinT *c, Grid *g, bool isFS); 

    virtual ~PFMSolver() {

        delete hypre;

        free_mem(soln_);
        free_mem(G1);
        free_mem(unk);

        if(grid->myid==0 && conf->verbose)
        printf("INFO: the memory allocate by PFMSolver has been released.\n");
    }

    virtual void evolve();         // evolve over a time slice 
    void newton_raphson(); // New-Raphson method iteration 
    
    void init();
    void init1d();
    void init2d();
    void init3d();

    inline void stencil() {
        if(ndim==3)      stencil_ac_3d_(grid->nxyz, lamdaxyz, &nguard, bcp, soln, &theta, &dtk, &beta_);
        else if(ndim==2) stencil_ac_2d_(grid->nxyz, lamdaxyz, &nguard, bcp, soln, &theta, &dtk, &beta_);
        else if(ndim==1) stencil_ac_1d_(grid->nxyz, lamdaxyz, &nguard, bcp, soln, &theta, &dtk, &beta_);
    }

    inline void rhs() {
        if(ndim==3)      rhs_ac_3d_(grid->nxyz, lamdaxyz, &nguard, b, soln, soln_, G1, &theta, &dtk, &beta_);
        else if(ndim==2) rhs_ac_2d_(grid->nxyz, lamdaxyz, &nguard, b, soln, soln_, G1, &theta, &dtk, &beta_);
        else if(ndim==1) rhs_ac_1d_(grid->nxyz, lamdaxyz, &nguard, b, soln, soln_, G1, &theta, &dtk, &beta_);
    }
    inline void rhs_g1() {
        if(ndim==3)      rhs_g1_ac_3d_(grid->nxyz, lamdaxyz, &nguard, soln_,  G1, &theta, &dtk, &beta_);
        else if(ndim==2) rhs_g1_ac_2d_(grid->nxyz, lamdaxyz, &nguard, soln_,  G1, &theta, &dtk, &beta_);
        else if(ndim==1) rhs_g1_ac_1d_(grid->nxyz, lamdaxyz, &nguard, soln_,  G1, &theta, &dtk, &beta_);
    }

    /*
     * NOTE : for space division, converge check must be performed in the whole geographical space 
     */
    void chk_eps(double *err); 

    inline void update() {
        if(ndim==3)      update_ac_3d_(grid->nxyz, &nguard, soln, unk);
        else if(ndim==2) update_ac_2d_(grid->nxyz, &nguard, soln, unk);
        else if(ndim==1) update_ac_1d_(grid->nxyz, &nguard, soln, unk);
    }
};

//problem specific parameter parsing
int pfm_inih(void* obj, const char* section, const char* name, const char* value);

#endif
