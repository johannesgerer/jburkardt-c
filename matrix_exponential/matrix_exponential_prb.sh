#!/bin/bash
#
gcc -c -g -I/$HOME/include matrix_exponential_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling matrix_exponential_prb.c."
  exit
fi
rm compiler.txt
#
gcc matrix_exponential_prb.o /$HOME/libc/$ARCH/matrix_exponential.o \
                             /$HOME/libc/$ARCH/test_matrix_exponential.o \
                             /$HOME/libc/$ARCH/c8lib.o \
                             /$HOME/libc/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading matrix_exponential_prb.o."
  exit
fi
#
rm matrix_exponential_prb.o
#
mv a.out matrix_exponential_prb
./matrix_exponential_prb > matrix_exponential_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running matrix_exponential_prb."
  exit
fi
rm matrix_exponential_prb
#
echo "Program output written to matrix_exponential_prb_output.txt"
