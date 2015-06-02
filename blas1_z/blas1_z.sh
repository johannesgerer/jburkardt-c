#!/bin/bash
#
cp blas1_z.h /$HOME/include
#
gcc -c -O2 -I$HOME/include blas1_z.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blas1_z.c."
  exit
fi
rm compiler.txt
#
mv blas1_z.o ~/libc/$ARCH/blas1_z.o
#
echo "Library installed as ~/libc/$ARCH/blas1_z.o"
