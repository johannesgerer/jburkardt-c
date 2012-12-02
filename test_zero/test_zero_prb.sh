#!/bin/bash
#
gcc -c -g -I/$HOME/include test_zero_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_zero_prb.c"
  exit
fi
rm compiler.txt
#
gcc test_zero_prb.o /$HOME/libc/$ARCH/test_zero.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_zero_prb.o."
  exit
fi
#
rm test_zero_prb.o
#
mv a.out test_zero_prb
./test_zero_prb > test_zero_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_zero_prb."
  exit
fi
rm test_zero_prb
#
echo "Program output written to test_zero_prb_output.txt"
