#ifndef _PinT_PFMPARAMS_H_
#define _PinT_PFMPARAMS_H_ 1

#include "PinT.h"
#include "Grid.h"

//problem specific parameter parsing
int pfm_inih(void* obj, const char* section, const char* name, const char* value);

/*
 *  NOTE :  
 *
 *  The steady condition is different for 1D/2D/3D due to the diffuse effect and reactive energy item. 
 *     for 1D : beta=0 .and.  any bc condition
 *     for 2/3D : still unknown 
 *
 *  So the final results are also different at 1D/2D/3D even if the program is absolutely correct.
 *  There is a simple method to verify the correctness of the program separately at 1D/2D/3D. 
 *  That is simulating the initial condition of 1D, 
 *  for 2D example, assigning the same value to the cell with the same x-coordinate: grid[i,...]=ini_val,
 *  that make all 1D vectors along y direction is the same (duplicated several 1D).
 *  Therefore the 2D result can be comparable with the 1D.
 *  
 *  BUT..., it is not enough to do above only, 
 *  the BC condition in Y direction will still have impact on the final result.
 *  The reason just is that diffuse effect will occur in Y direction if the neighbour points have different value.
 *  Thus, the 2D cannot completely simulate the 1D, their final result is no longer similar.
 *  You can use the customized BC function provided by the framework to implement any conditions.  
 *
 *  As a summary, numerical analysis may become much more complicated when being applied to real world problems. 
 *  That is a really hard work !!!
 *
 * *********************************
 *
 *  WARN : 
 *
 *  Nonlinear easily lead to unconvergence of linear solver. 
 *  Though implict algorithms usually don't require a smalll timestep,    
 *  if the timestep is big enough, the nonlinear item will also cause unconvergence or unphysical result.    
 *  The reaction item of the following Allen-Cahn equation, has a max value nearly 0.05 when u in [0,1] region if beta=0 and k=1  
 *  If making sure convergence, the value of dt*k*0.05 is better less than ONE, so the dtk=dt*k is should be less 20. 
 *  According our experiments, if the dtk parameter is bigger than 1.6, newton_raphson method is easily broken. 
 * 
 *  d{u}/d{t} = D*(Laplace operator){u} - k*{u}*({u}-1.0)*({u}-0.5+beta)
 *    {u} : unknown variables
 *    {t} : time
 *     D  : diffusion coefficient
 *     k  : interfacial width related param, it is more bigger, the interfacial width is more thiner 
 *
 * The convergence of time direction will become more difficult from 1D to 3D if there isn't a steady state of the equation.
 * It is hard to find the proper paramters that can quickly lead to a heat equilibrium at 3D than 1D/2D. 
 * Bigger timestep will also lead Crank-Nicolson method to a more divergent initial result before the first time iteration, 
 * that will absolutely cause the following calculations much more divergent.  
 *
 * Let lamda = D*dt/{dx^2}, dx=dy=dz
 * The CFL condition is lamda<1/2  at 1D, lamda<1/4 at 2D, and lamda<1/6 at 3D 
 * Therefore, even the number of space division is same in 1D and 3D, the max timestep in 3D becomes 1/3 of that in 1D.
 * max timestep = 1/6*dx^2
 * We further reduced the theoretical max timestep of explict solver by a factor of 0.311 for both convergence and testing convenience. 
 * Thus in a 144^3 unit cube, dt=2.5e-6 and the total timesteps=20k (serial fine solver) 
 *
 */
struct PFMParams {

    double k;     // interfacial width related
    double d;     // diffuse coefficient 
    double beta;  // potential engrgy parameter
    double beta_;       // 0.5-beta

    double theta = 0.5;  //Crank-Nicolson; 0: Ex.Euler; 1: Im.Euler

    double dtk;         // dt*k
    double lamda_x;     // d*dt/(dx**2)
    double lamda_y = 0.0;     
    double lamda_z = 0.0;    
    double lamdaxyz[3];

    int newton_itmax = 10;
    
    int ndim; 

    int fsolver = 0;
    int csolver = 0;
    
    void init(PinT* conf, Grid *g){
        conf->init_module(this, pfm_inih);
        beta_ = 0.5 - beta;
        ndim = conf->ndim;
    }

    void init(PinT* conf, Grid *g, bool isFS){
        init(conf, g);

        if(isFS) {
            lamda_x = d*conf->f_dt/(g->dx*g->dx);   
            if(ndim>=2) lamda_y = d * conf->f_dt / (g->dy*g->dy);
            if(ndim>=3) lamda_z = d * conf->f_dt / (g->dz*g->dz);
            dtk = conf->f_dt*k;
        }else {
            lamda_x = d*conf->c_dt/(g->dx*g->dx);   
            if(ndim>=2) lamda_y = d * conf->c_dt / (g->dy*g->dy);
            if(ndim>=3) lamda_z = d * conf->c_dt / (g->dz*g->dz);
            dtk = conf->c_dt*k; 
        }

        lamdaxyz[0] = lamda_x;
        lamdaxyz[1] = lamda_y;
        lamdaxyz[2] = lamda_z;
    }

    void print() {
        printf("PFM init parameter : \n");  
        printf("        D    = %f \n", d);
        printf("        k    = %f \n", k);
        printf("        beta = %f \n", beta);
        printf("  newton_iter= %d \n", newton_itmax);
        printf("        theta= %f , (0.0: Ex.Euler; 0.5: Crank-Nicolson; 1.0: Im.Euler.)\n\n", theta);
        printf("  Fine Solver= %s \n", getSolvername(fsolver));
        printf("  Coar Solver= %s \n", getSolvername(csolver));
    
    }
    void printLamda(const char* flag) {
        printf("%s : dtk=%f, lamda_x=%f, lamday=%f, lamdaz=%f \n", flag, dtk, lamda_x, lamda_y, lamda_z );
    }

    const char* CN  = "CN";
    const char* BD  = "BD4";
    const char* EU  = "Euler";
    const int ID_CN = 0;
    const int ID_BD = 1;
    const int ID_EU = 2;

    const char* getSolvername(int s) {
        if(s==ID_CN) return CN;
        if(s==ID_BD) return BD;
        if(s==ID_EU) return EU;
    }
};

#endif

