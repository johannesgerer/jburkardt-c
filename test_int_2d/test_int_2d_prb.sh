#!/bin/bash
#
gcc -c -g -I/$HOME/include test_int_2d_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_int_2d_prb.c"
  exit
fi
rm compiler.txt
#
gcc test_int_2d_prb.o /$HOME/libc/$ARCH/test_int_2d.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_int_2d_prb.o."
  exit
fi
#
rm test_int_2d_prb.o
#
mv a.out test_int_2d_prb
./test_int_2d_prb > test_int_2d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_int_2d_prb."
  exit
fi
rm test_int_2d_prb
#
echo "Program output written to test_int_2d_prb_output.txt"
