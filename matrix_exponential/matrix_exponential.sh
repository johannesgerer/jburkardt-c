#!/bin/bash
#
cp matrix_exponential.h /$HOME/include
#
gcc -c -g -I/$HOME/include matrix_exponential.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling matrix_exponential.c."
  exit
fi
rm compiler.txt
#
mv matrix_exponential.o ~/libc/$ARCH/matrix_exponential.o
#
echo "Library installed as ~/libc/$ARCH/matrix_exponential.o"
