#!/bin/bash
#
gcc -c -g -I/$HOME/include qr_solve_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling qr_solve_prb.c"
  exit
fi
rm compiler.txt
#
gcc qr_solve_prb.o /$HOME/libc/$ARCH/qr_solve.o /$HOME/libc/$ARCH/test_ls.o  /$HOME/libc/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading qr_solve_prb.o."
  exit
fi
#
rm qr_solve_prb.o
#
mv a.out qr_solve_prb
./qr_solve_prb > qr_solve_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running qr_solve_prb."
  exit
fi
rm qr_solve_prb
#
echo "Program output written to qr_solve_prb_output.txt"
