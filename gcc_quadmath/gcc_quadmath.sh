#!/bin/bash
#
gcc -c gcc_quadmath.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling gcc_quadmath.c."
  exit
fi
rm compiler.txt
#
gcc gcc_quadmath.o
if [ $? -ne 0 ]; then
  echo "Errors linking gcc_quadmath.o."
  exit
fi
#
rm gcc_quadmath.o
#
mv a.out gcc_quadmath
./gcc_quadmath > gcc_quadmath_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running gcc_quadmath."
  exit
fi
rm gcc_quadmath
#
echo "Program output written to gcc_quadmath_output.txt"
