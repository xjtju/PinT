#ifndef PinT_DRIVER_H_
#define PinT_DRIVER_H_ 1

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <mpi.h>

#include "ini.h"
#include "PinT.h"

/**
 * the driver of the PinT process 
 */

class Driver {

    int myid = 0;
 public:

    void init(char* ini_file, PinT* conf);

    void evolve();

    void finalize();

    void INFO (const char* fmt, ...);
    void WARN (const char* fmt, ...);
    void Abort(const char* fmt, ...);

};

static int handler(void* pint, const char* section, const char* name, const char* value);

#endif
