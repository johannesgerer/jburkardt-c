#!/bin/bash
#
gcc -c -g -I/$HOME/include qwgw_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling qwgw_prb.c"
  exit
fi
rm compiler.txt
#
gcc qwgw_prb.o /$HOME/libc/$ARCH/qwgw.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading qwgw_prb.o."
  exit
fi
#
rm qwgw_prb.o
#
mv a.out qwgw_prb
./qwgw_prb > qwgw_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running qwgw_prb."
  exit
fi
rm qwgw_prb
#
echo "Program output written to qwgw_prb_output.txt"
