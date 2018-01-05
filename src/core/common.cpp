#include "common.h"

double* alloc_mem(size_t size){

    double *u;  
    if( !(u = new double[size]) ) {
        printf("Failed to allocate memory %d [MB]\n", size/MB );
        exit(1); 
    }
    for(int i=0; i<size; i++) u[i] = 0.0;

    return u;
}

void free_mem(double *u){
    delete u;
}
