#!/bin/bash
#
#  Compile the program with GCC.
#
/usr/local/bin/gcc -fopenmp helmholtz.c -lm
#
#  Compile the program with ICC.
#
#icc -openmp -parallel helmholtz.c -lm
#
mv a.out helmholtz
#
#  Run with 1, 2, and 4 threads.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./helmholtz > helmholtz_local_output.txt
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./helmholtz >> helmholtz_local_output.txt
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./helmholtz >> helmholtz_local_output.txt
#
#  Discard the executable file.
#
rm helmholtz
#
echo "Program output written to helmholtz_local_output.txt"
