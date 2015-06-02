#!/bin/bash
#
cp laplacian.h /$HOME/include
#
gcc -c -g -I/$HOME/include laplacian.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling laplacian.c"
  exit
fi
rm compiler.txt
#
mv laplacian.o ~/libc/$ARCH/laplacian.o
#
echo "Library installed as ~/libc/$ARCH/laplacian.o"
