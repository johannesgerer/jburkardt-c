#!/bin/bash
#
#  Compile the program with GCC.
#
/usr/local/bin/gcc -fopenmp quad2d_openmp.c -lm
#
mv a.out quad2d
#
#  Run with 1, 2, and 4 threads.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./quad2d > quad2d_local_gcc_output.txt
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./quad2d >> quad2d_local_gcc_output.txt
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./quad2d >> quad2d_local_gcc_output.txt
#
#  Discard the executable file.
#
rm quad2d
#
echo "Program output written to quad2d_local_gcc_output.txt"
