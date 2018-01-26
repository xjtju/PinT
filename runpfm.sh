export OMP_NUM_THREADS=1
mpirun -host localhost -np 1 ./pfm_alpha.exe pfm.ini pfm.bd4
