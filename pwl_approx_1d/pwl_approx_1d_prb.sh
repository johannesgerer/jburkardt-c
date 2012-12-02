#!/bin/bash
#
gcc -c -g -I/$HOME/include pwl_approx_1d_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pwl_approx_1d_prb.c"
  exit
fi
rm compiler.txt
#
gcc pwl_approx_1d_prb.o /$HOME/libc/$ARCH/pwl_approx_1d.o \
                        /$HOME/libc/$ARCH/test_interp_1d.o \
                        /$HOME/libc/$ARCH/qr_solve.o \
                        /$HOME/libc/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pwl_approx_1d_prb.o."
  exit
fi
#
rm pwl_approx_1d_prb.o
#
mv a.out pwl_approx_1d_prb
./pwl_approx_1d_prb > pwl_approx_1d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pwl_approx_1d_prb."
  exit
fi
rm pwl_approx_1d_prb
#
echo "Program output written to pwl_approx_1d_prb_output.txt"
