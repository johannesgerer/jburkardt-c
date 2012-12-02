#!/bin/bash
#
gcc -c -g -I/$HOME/include lagrange_approx_1d_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lagrange_approx_1d_prb.c"
  exit
fi
rm compiler.txt
#
gcc lagrange_approx_1d_prb.o /$HOME/libc/$ARCH/lagrange_approx_1d.o \
                             /$HOME/libc/$ARCH/test_interp_1d.o \
                             /$HOME/libc/$ARCH/qr_solve.o \
                             /$HOME/libc/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lagrange_approx_1d_prb.o."
  exit
fi
#
rm lagrange_approx_1d_prb.o
#
mv a.out lagrange_approx_1d_prb
./lagrange_approx_1d_prb > lagrange_approx_1d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lagrange_approx_1d_prb."
  exit
fi
rm lagrange_approx_1d_prb
#
echo "Program output written to lagrange_approx_1d_prb_output.txt"
