#!/bin/bash
#
gcc -c -g -I/$HOME/include qwv_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling qwv_prb.c"
  exit
fi
rm compiler.txt
#
gcc qwv_prb.o /$HOME/libc/$ARCH/qwv.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading qwv_prb.o."
  exit
fi
#
rm qwv_prb.o
#
mv a.out qwv_prb
./qwv_prb > qwv_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running qwv_prb."
  exit
fi
rm qwv_prb
#
echo "Program output written to qwv_prb_output.txt"
