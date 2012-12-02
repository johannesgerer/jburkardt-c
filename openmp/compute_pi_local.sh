#!/bin/bash
#
#  Compile the program with GCC.
#
/usr/local/bin/gcc -fopenmp compute_pi.c -lm
#
#  Compile the program with ICC.
#
#icc -openmp -parallel compute_pi.c -lm
#
mv a.out compute_pi
#
#  Run with 1, 2, and 4 threads.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./compute_pi > compute_pi_local_output.txt
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./compute_pi >> compute_pi_local_output.txt
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./compute_pi >> compute_pi_local_output.txt
#
#  Discard the executable file.
#
rm compute_pi
#
echo "Program output written to compute_pi_local_output.txt"
