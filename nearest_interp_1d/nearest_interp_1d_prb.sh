#!/bin/bash
#
gcc -c -g -I/$HOME/include nearest_interp_1d_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nearest_interp_1d_prb.c"
  exit
fi
rm compiler.txt
#
gcc nearest_interp_1d_prb.o /$HOME/libc/$ARCH/nearest_interp_1d.o \
                            /$HOME/libc/$ARCH/test_interp.o \
                            /$HOME/libc/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading nearest_interp_1d_prb.o"
  exit
fi
#
rm nearest_interp_1d_prb.o
#
mv a.out nearest_interp_1d_prb
./nearest_interp_1d_prb > nearest_interp_1d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running nearest_interp_1d_prb."
  exit
fi
rm nearest_interp_1d_prb
#
echo "Program output written to nearest_interp_1d_prb_output.txt"
