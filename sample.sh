mpicc pfc.c -lfftw3_mpi -lfftw3 -lm -Ofast -Wall -o pfc
mpirun -np 8 pfc step1
mpirun -np 8 pfc step2
java -jar plotter.jar step1-t:# step1-t:# 1024 1024 0 1000 10000
java -jar plotter.jar step2-t:# step2-t:# 1024 1024 0 1000 10000
