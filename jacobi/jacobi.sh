#!/bin/bash
#
cp jacobi.h /$HOME/include
#
gcc -c -g -I /$HOME/include jacobi.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling jacobi.c"
  exit
fi
rm compiler.txt
#
mv jacobi.o ~/libc/$ARCH/jacobi.o
#
echo "Library installed as ~/libc/$ARCH/jacobi.o"
