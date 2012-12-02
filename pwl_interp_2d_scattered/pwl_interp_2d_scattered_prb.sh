#!/bin/bash
#
gcc -c -g -I/$HOME/include pwl_interp_2d_scattered_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pwl_interp_2d_scattered_prb.c"
  exit
fi
rm compiler.txt
#
gcc pwl_interp_2d_scattered_prb.o /$HOME/libc/$ARCH/pwl_interp_2d_scattered.o \
                                  /$HOME/libc/$ARCH/test_interp_2d.o \
                                  /$HOME/libc/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pwl_interp_2d_scattered_prb.o."
  exit
fi
#
rm pwl_interp_2d_scattered_prb.o
#
mv a.out pwl_interp_2d_scattered_prb
./pwl_interp_2d_scattered_prb > pwl_interp_2d_scattered_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pwl_interp_2d_scattered_prb."
  exit
fi
rm pwl_interp_2d_scattered_prb
#
echo "Program output written to pwl_interp_2d_scattered_prb_output.txt"
