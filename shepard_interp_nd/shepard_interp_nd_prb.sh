#!/bin/bash
#
gcc -c -g -I/$HOME/include shepard_interp_nd_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling shepard_interp_nd_prb.c"
  exit
fi
rm compiler.txt
#
gcc shepard_interp_nd_prb.o /$HOME/libc/$ARCH/shepard_interp_nd.o \
                            /$HOME/libc/$ARCH/test_interp_nd.o \
                            /$HOME/libc/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading shepard_interp_nd_prb.o."
  exit
fi
#
rm shepard_interp_nd_prb.o
#
mv a.out shepard_interp_nd_prb
./shepard_interp_nd_prb > shepard_interp_nd_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running shepard_interp_nd_prb."
  exit
fi
rm shepard_interp_nd_prb
#
echo "Program output written to shepard_interp_nd_prb_output.txt"
