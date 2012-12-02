#!/bin/bash
#
gcc -c -g -I/$HOME/include test_min_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_min_prb.c"
  exit
fi
rm compiler.txt
#
gcc test_min_prb.o /$HOME/libc/$ARCH/test_min.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_min_prb.o."
  exit
fi
#
rm test_min_prb.o
#
mv a.out test_min_prb
./test_min_prb > test_min_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_min_prb."
  exit
fi
rm test_min_prb
#
echo "Program output written to test_min_prb_output.txt"
