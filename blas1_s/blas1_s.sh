#!/bin/bash
#
cp blas1_s.h /$HOME/include
#
gcc -c -O2 -I$HOME/include blas1_s.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blas1_s.c."
  exit
fi
rm compiler.txt
#
mv blas1_s.o ~/libc/$ARCH/blas1_s.o
#
echo "Library installed as ~/libc/$ARCH/blas1_s.o"
