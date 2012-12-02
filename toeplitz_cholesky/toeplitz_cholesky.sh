#!/bin/bash
#
cp toeplitz_cholesky.h /$HOME/include
#
gcc -c -g -I/$HOME/include toeplitz_cholesky.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toeplitz_cholesky.c"
  exit
fi
rm compiler.txt
#
mv toeplitz_cholesky.o ~/libc/$ARCH/toeplitz_cholesky.o
#
echo "Library installed as ~/libc/$ARCH/toeplitz_cholesky.o"
