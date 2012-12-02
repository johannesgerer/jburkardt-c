#!/bin/bash
#
gcc -c -g -I/$HOME/include barycentric_interp_1d_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling barycentric_interp_1d_prb.c"
  exit
fi
rm compiler.txt
#
gcc barycentric_interp_1d_prb.o /$HOME/libc/$ARCH/barycentric_interp_1d.o \
                                /$HOME/libc/$ARCH/test_interp_1d.o \
                                /$HOME/libc/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading barycentric_interp_1d_prb.o."
  exit
fi
#
rm barycentric_interp_1d_prb.o
#
mv a.out barycentric_interp_1d_prb
./barycentric_interp_1d_prb > barycentric_interp_1d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running barycentric_interp_1d_prb."
  exit
fi
rm barycentric_interp_1d_prb
#
echo "Program output written to barycentric_interp_1d_prb_output.txt"
