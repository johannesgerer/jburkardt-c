#!/bin/bash
#
#  Compile the program.
#
gcc ../hello_openmp/hello_openmp.c -I. -L$HOME/libc/$ARCH -lopenmp_stubs
mv a.out hello
#
#  Run the program.
#
./hello > hello_output.txt
rm hello
#
echo "Normal end of execution."
