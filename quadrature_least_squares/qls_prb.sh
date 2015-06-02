#!/bin/bash
#
gcc -c -I/$HOME/include qls_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling qls_prb.c"
  exit
fi
rm compiler.txt
#
gcc qls_prb.o /$HOME/libc/$ARCH/qls.o /$HOME/libc/$ARCH/qr_solve.o /$HOME/libc/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading qls_prb.o."
  exit
fi
#
rm qls_prb.o
#
mv a.out qls_prb
./qls_prb > qls_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running qls_prb."
  exit
fi
rm qls_prb
#
echo "Program output written to qls_prb_output.txt"
