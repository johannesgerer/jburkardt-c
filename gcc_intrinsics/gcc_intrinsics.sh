#!/bin/bash
#
gcc -c gcc_intrinsics.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling gcc_intrinsics.c."
  exit
fi
rm compiler.txt
#
gcc gcc_intrinsics.o
if [ $? -ne 0 ]; then
  echo "Errors linking gcc_intrinsics.o."
  exit
fi
#
rm gcc_intrinsics.o
#
mv a.out gcc_intrinsics
./gcc_intrinsics > gcc_intrinsics_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running gcc_intrinsics."
  exit
fi
rm gcc_intrinsics
#
echo "Program output written to gcc_intrinsics_output.txt"
