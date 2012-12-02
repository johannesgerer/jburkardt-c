#!/bin/bash
#
#  Compile the program.
#
gcc ../openmp/compute_pi.c -I. -L$HOME/libc/$ARCH -lopenmp_stubs
mv a.out compute_pi
#
#  Run the program.
#
./compute_pi > compute_pi_output.txt
rm compute_pi
#
echo "Normal end of execution."
