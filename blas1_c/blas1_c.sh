#!/bin/bash
#
cp blas1_c.h /$HOME/include
#
gcc -c -O2 -I/$HOME/include blas1_c.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blas1_c.c."
  exit
fi
rm compiler.txt
#
mv blas1_c.o ~/libc/$ARCH/blas1_c.o
#
echo "Library installed as ~/libc/$ARCH/blas1_c.o"
