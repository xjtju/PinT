export OMP_NUM_THREADS=1
mpirun -host localhost -np 1 ./pfm_alpha.exe src/pfm/pfm.ini pfm.bd4
