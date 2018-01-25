export OMP_NUM_THREADS=1
mpirun -host localhost -np 4 ./pfm_alpha.exe pfm.ini pfm.32.4
