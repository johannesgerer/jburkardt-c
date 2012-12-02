#!/bin/bash
#
#  Compile the program with GCC.
#
gcc -fopenmp prime_openmp.c
#
#  Compile the program with ICC.
#
#icc -openmp -parallel prime_openmp.c
#
mv a.out prime
#
#  Run with 1, 2, 4 and 8 threads.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./prime > prime_local_output.txt
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./prime >> prime_local_output.txt
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./prime >> prime_local_output.txt
#
echo "Run with 8 threads."
export OMP_NUM_THREADS=8
./prime >> prime_local_output.txt
#
#  Discard the executable file.
#
rm prime
#
echo "Program output written to prime_local_output.txt"
