#!/bin/bash
#
gcc -c quadrature.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling quadrature.c"
  exit
fi
rm compiler.txt
#
gcc quadrature.o $HOME/libc/$ARCH/mpi_stubs.o
if [ $? -ne 0 ]; then
  echo "Errors while linking and loading quadrature.o"
  exit
fi
rm quadrature.o
#
mv a.out quadrature
./quadrature > quadrature_output.txt
if [ $? -ne 0 ]; then
  echo "Errors while running quadrature"
  exit
fi
rm quadrature
#
echo "Program output written to quadrature_output.txt"
