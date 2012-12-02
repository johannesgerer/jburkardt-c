#!/bin/bash
#
/usr/local/bin/gcc -fopenmp timer_omp_get_wtime.c
mv a.out timer_omp_get_wtime
#
#  Run the program.
#
./timer_omp_get_wtime > timer_omp_get_wtime_output.txt
rm timer_omp_get_wtime
#
echo "Program output written to timer_omp_get_wtime_output.txt"
