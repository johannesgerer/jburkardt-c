#!/bin/bash
#
cp blas1_d.h /$HOME/include
#
gcc -c -O2 -I$HOME/include blas1_d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blas1_d.c."
  exit
fi
rm compiler.txt
#
mv blas1_d.o ~/libc/$ARCH/blas1_d.o
#
echo "Library installed as ~/libc/$ARCH/blas1_d.o"
