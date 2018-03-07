
struct BiCGParam {


    int itmax = 20; 
    double eps = 1.0e-6;
    int force_abort = 0; 
    
    bool isPrecond = false; // in the current version, no preconditioner is BETTER according to experiments.

    int sor_itmax = 4;
    double sor_omega = 1.5;
    double sor_checkCnvg = false;
    
    void init(PinT* conf){
        conf->init_module(this, BiCG_inih);
    }
}
