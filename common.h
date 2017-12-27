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

double* alloc_mem(size_t size);
void clear_mem(double *x, size_t size);
void free_mem(double *u);

#endif
