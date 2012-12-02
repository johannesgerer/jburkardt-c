#!/bin/bash
#
gcc -c -g -I/$HOME/include fn_prb2.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fn_prb2.c"
  exit
fi
rm compiler.txt
#
gcc fn_prb2.o /$HOME/libc/$ARCH/fn.o /$HOME/libc/$ARCH/test_values.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fn_prb2.o."
  exit
fi
#
rm fn_prb2.o
#
mv a.out fn_prb2
./fn_prb2 > fn_prb2_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fn_prb2."
  exit
fi
rm fn_prb2
#
echo "Program output written to fn_prb2_output.txt"
