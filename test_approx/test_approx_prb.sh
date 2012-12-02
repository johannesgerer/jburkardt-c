#!/bin/bash
#
gcc -c -g -I/$HOME/include test_approx_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_approx_prb.c"
  exit
fi
rm compiler.txt
#
gcc test_approx_prb.o /$HOME/libc/$ARCH/test_approx.o /$HOME/libc/$ARCH/spline.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_approx_prb.o."
  exit
fi
#
rm test_approx_prb.o
#
mv a.out test_approx_prb
./test_approx_prb > test_approx_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_approx_prb."
  exit
fi
rm  test_approx_prb
#
echo "Program output written to test_approx_prb_output.txt"
