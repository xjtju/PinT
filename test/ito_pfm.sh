#!/bin/bash

#PJM -L "rscunit=ito-a"
#PJM -L "rscgrp=ito-l"
#PJM -L "vnode=1600"
#PJM -L "vnode-core=1"
#PJM -L "elapse=04:00:00"
#PJM -e "ito.perf"
#PJM -o "ito.info"
#PJM -X

# OpenMP parallel thread number
#export OMP_NUM_THREADS=18

# automatically parallel thread number
#export PARALLEL=18 
module load hdf5/1.10.1-non-parallel-fj 
# argv[1] : ini file
# argv[2] : perfix of output file name, not necessary   

mpiexec -n 1600 ../pfm.exe ./fpfm.ini
