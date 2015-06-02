#!/bin/bash
#
gcc -c test_mat_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling test_mat_prb.c."
  exit
fi
#
gcc test_mat_prb.o /$HOME/libc/$ARCH/test_mat.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_mat_prb.o."
  exit
fi
#
rm test_mat_prb.o
#
mv a.out test_mat_prb
./test_mat_prb > test_mat_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_mat_prb."
  exit
fi
rm test_mat_prb
#
echo "Program output written to test_mat_prb_output.txt"
