#!/bin/bash
#
#  Compile the program.
#
gcc ../openmp/helmholtz.c -I. -L$HOME/libc/$ARCH -lopenmp_stubs
mv a.out helmholtz
#
#  Run the program.
#
./helmholtz > helmholtz_output.txt
rm helmholtz
#
echo "Normal end of execution."
