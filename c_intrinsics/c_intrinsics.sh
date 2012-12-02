#!/bin/bash
#
gcc -c c_intrinsics.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling c_intrinsics.c."
  exit
fi
rm compiler.txt
#
gcc c_intrinsics.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking c_intrinsics.o."
  exit
fi
#
rm c_intrinsics.o
#
mv a.out c_intrinsics
./c_intrinsics > c_intrinsics_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running c_intrinsics."
  exit
fi
rm c_intrinsics
#
echo "Program output written to c_intrinsics_output.txt"
