#!/bin/bash
#
cp super_blas.h /$HOME/include
#
gcc -c -g super_blas.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling super_blas.c."
  exit
fi
rm compiler.txt
#
mv super_blas.o ~/libc/$ARCH/super_blas.o
#
echo "Library installed as ~/libc/$ARCH/super_blas.o"
