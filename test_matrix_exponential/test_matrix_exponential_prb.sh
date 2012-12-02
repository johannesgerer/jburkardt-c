#!/bin/bash
#
gcc -c -g -I/$HOME/include test_matrix_exponential_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_matrix_exponential_prb.c."
  exit
fi
rm compiler.txt
#
gcc test_matrix_exponential_prb.o /$HOME/libc/$ARCH/test_matrix_exponential.o \
                                  /$HOME/libc/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_matrix_exponential_prb.o."
  exit
fi
#
rm test_matrix_exponential_prb.o
#
mv a.out test_matrix_exponential_prb
./test_matrix_exponential_prb > test_matrix_exponential_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_matrix_exponential_prb."
  exit
fi
rm test_matrix_exponential_prb
#
echo "Program output written to test_matrix_exponential_prb_output.txt"
