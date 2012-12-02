#!/bin/bash
#
gcc -c c_operators.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling c_operators.c."
  exit
fi
rm compiler.txt
#
gcc c_operators.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking c_operators.o."
  exit
fi
#
rm c_operators.o
#
mv a.out c_operators
./c_operators > c_operators_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running c_operators."
  exit
fi
rm c_operators
#
echo "Program output written to c_operators_output.txt"
