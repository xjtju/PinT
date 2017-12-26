/**
 * PinT -- Hello world 
 * 
 * Xiao Jian, 2017.12.22
 */
#ifndef PinT_COMMON_H
#define PinT_COMMON_H 1

#include <math.h>
#include <assert.h>
#include <stdio.h>

const double PI =  M_PI;
const long MB = 1024*1024;

//the size of time-space domain
const double Tspan = 4.0;
const double Xspan = 1.0; 
const long NT = 10000;
const int  NX = 100;
const int NGuard = 1; //guard cells

// the steps of time-space domain
const double DT = Tspan/NT; 
const double DX = Xspan / (NX+2*NGuard-1);

//the number of time slices, parallel processes along time domain 
const int Tslices = 10;
const int Tsteps_oneslice = NT/Tslices; 

const double EPS = 1.0e-6;

double* alloc_mem(size_t size);
void clear_mem(double *x, size_t size);
void free_mem(double *u);

#endif
