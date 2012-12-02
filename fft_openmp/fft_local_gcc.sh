#!/bin/bash
#
#  Compile the program with GCC.
#
/usr/local/bin/gcc -fopenmp fft_openmp.c -lm
#
mv a.out fft
#
#  Run with 1, 2, and 4 threads.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./fft > fft_local_gcc_output.txt
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./fft >> fft_local_gcc_output.txt
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./fft >> fft_local_gcc_output.txt
#
#  Discard the executable file.
#
rm fft
#
echo "Program output written to fft_local_gcc_output.txt"
