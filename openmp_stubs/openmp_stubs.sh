#!/bin/bash
#
cp omp.h /$HOME/include
#
gcc -c -g openmp_stubs.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling openmp_stubs.c."
  exit
fi
rm compiler.txt
#
ar qc libopenmp_stubs.a openmp_stubs.o
mv libopenmp_stubs.a ~/libc/$ARCH/libopenmp_stubs.a
rm openmp_stubs.o
#
echo "Library installed as ~/libc/$ARCH/libopenmp_stubs.a"
