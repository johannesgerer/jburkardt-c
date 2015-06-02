#!/bin/bash
#
cp blas3.h /$HOME/include
#
gcc -c -O2 -I$HOME/include blas3.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blas3.c."
  exit
fi
rm compiler.txt
#
mv blas3.o ~/libc/$ARCH/blas3.o
#
echo "Library installed as ~/libc/$ARCH/blas3.o"
