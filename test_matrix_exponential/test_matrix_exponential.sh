#!/bin/bash
#
cp test_matrix_exponential.h /$HOME/include
#
gcc -c -g -I /$HOME/include test_matrix_exponential.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_matrix_exponential.c."
  exit
fi
rm compiler.txt
#
mv test_matrix_exponential.o ~/libc/$ARCH/test_matrix_exponential.o
#
echo "Library installed as ~/libc/$ARCH/test_matrix_exponential.o"
