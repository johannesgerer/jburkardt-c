#!/bin/bash
#
gcc -c -g -I/$HOME/include vandermonde_approx_2d_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling vandermonde_approx_2d_prb.c"
  exit
fi
rm compiler.txt
#
gcc vandermonde_approx_2d_prb.o /$HOME/libc/$ARCH/vandermonde_approx_2d.o \
                                /$HOME/libc/$ARCH/test_interp_2d.o \
                                /$HOME/libc/$ARCH/qr_solve.o \
                                /$HOME/libc/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading vandermonde_approx_2d_prb.o."
  exit
fi
#
rm vandermonde_approx_2d_prb.o
#
mv a.out vandermonde_approx_2d_prb
./vandermonde_approx_2d_prb > vandermonde_approx_2d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running vandermonde_approx_2d_prb."
  exit
fi
rm vandermonde_approx_2d_prb
#
echo "Program output written to vandermonde_approx_2d_prb_output.txt"
