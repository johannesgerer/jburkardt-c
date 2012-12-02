#!/bin/bash
#
#  Compile the program with ICC.
#
icc -openmp -parallel schedule_openmp.c
#
mv a.out schedule
#
#  Run with 2 threads.
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./schedule > schedule_local_icc_output.txt
#
#  Discard the executable file.
#
rm schedule
#
echo "Program output written to schedule_local_icc_output.txt"
