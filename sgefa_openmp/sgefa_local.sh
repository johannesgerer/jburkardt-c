#!/bin/bash
#
#  Compile the program with GCC.
#
/usr/local/bin/gcc -fopenmp sgefa_openmp.c
#
#  Compile the program with ICC.
#
#icc -openmp -parallel sgefa_openmp.c
#
mv a.out sgefa
#
#  Run with 1, 2, 4 and 8 threads.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./sgefa > sgefa_local_output.txt
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./sgefa >> sgefa_local_output.txt
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./sgefa >> sgefa_local_output.txt
#
echo "Run with 8 threads."
export OMP_NUM_THREADS=8
./sgefa >> sgefa_local_output.txt
#
#  Discard the executable file.
#
rm sgefa
#
echo "Program output written to sgefa_local_output.txt"
