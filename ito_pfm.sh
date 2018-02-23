#!/bin/bash

#PJM -L "rscunit=ito-a"
#PJM -L "rscgrp=ito-ss-dbg"
#PJM -L "vnode=4"
#PJM -L "vnode-core=1"
#PJM -L "elapse=15:00"
#PJM -j
#PJM -X

# OpenMP parallel thread number
export OMP_NUM_THREADS=1

# automatically parallel thread number
#export PARALLEL=36 
module load hdf5/1.10.1-non-parallel-fj 
# argv[1] : ini file
# argv[2] : perfix of output file name  

mpiexec -n 4 ./pint_pfm.exe ./pfm.ini pfm.cn.0 
#mpiexec -n 36 ./pint_pfm.exe src/pfm/pfm.ini pfm.cn.1 
