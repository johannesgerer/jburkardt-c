#!/bin/bash
#
gcc -c -g -I/$HOME/include test_interp_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_interp_prb.c"
  exit
fi
rm compiler.txt
#
gcc test_interp_prb.o /$HOME/libc/$ARCH/test_interp.o /$HOME/libc/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_interp_prb.o."
  exit
fi
#
rm test_interp_prb.o
#
mv a.out test_interp_prb
./test_interp_prb > test_interp_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_interp_prb."
  exit
fi
rm test_interp_prb
#
echo "Program output written to test_interp_prb_output.txt"
