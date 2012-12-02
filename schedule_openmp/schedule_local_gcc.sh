#!/bin/bash
#
#  Compile the program with GCC.
#
/usr/local/bin/gcc -fopenmp schedule_openmp.c
#
mv a.out schedule
#
#  Run with 2 threads.
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./schedule > schedule_local_gcc_output.txt
#
#  Discard the executable file.
#
rm schedule
#
echo "Program output written to schedule_local_gcc_output.txt"
