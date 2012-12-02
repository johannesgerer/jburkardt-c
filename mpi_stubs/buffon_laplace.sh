#!/bin/bash
#
gcc -c buffon_laplace.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling buffon_laplace.c"
  exit
fi
rm compiler.txt
#
gcc buffon_laplace.o $HOME/libc/$ARCH/mpi_stubs.o -lm
if [ $? -ne 0 ]; then
  echo "Errors while linking and loading buffon_laplace.o"
  exit
fi
rm buffon_laplace.o
#
mv a.out buffon_laplace
./buffon_laplace > buffon_laplace_output.txt
if [ $? -ne 0 ]; then
  echo "Errors while running buffon_laplace"
  exit
fi
rm buffon_laplace
#
echo "Program output written to buffon_laplace_output.txt"
