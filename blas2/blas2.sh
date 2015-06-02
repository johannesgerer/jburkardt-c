#!/bin/bash
#
cp blas2.h /$HOME/include
#
gcc -c -O2 -I$HOME/include blas2.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blas2.c."
  exit
fi
rm compiler.txt
#
mv blas2.o ~/libc/$ARCH/blas2.o
#
echo "Library installed as ~/libc/$ARCH/blas2.o"
