#include "Model.h"
/**
 * the 1D Heat Diffusion with reflect boundary condition
 */
struct HeatModel : public Model {

public: 

    HeatModel():Model() { }

    // nguard = 1
    void bc(){
        x[0] = x[1];
        x[nx+nguard] = x[nx];
    }

    //Ut = kUxx, k = .061644
    int init_x(){
        for(int i = nguard; i<nx+nguard ; i++){
           x[i] = cos(2*i*DX); 
        }
        bc();
        return 0;
    }
};

