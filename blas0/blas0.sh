#!/bin/bash
#
cp blas0.h /$HOME/include
#
gcc -c -O2 blas0.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blas0.c."
  exit
fi
rm compiler.txt
#
mv blas0.o ~/libc/$ARCH/blas0.o
#
echo "Library installed as ~/libc/$ARCH/blas0.o"
