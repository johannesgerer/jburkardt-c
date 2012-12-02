#!/bin/bash
#
gcc -c -g -I/$HOME/include chebyshev_interp_1d_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling chebyshev_interp_1d_prb.c"
  exit
fi
rm compiler.txt
#
gcc chebyshev_interp_1d_prb.o /$HOME/libc/$ARCH/chebyshev_interp_1d.o \
                              /$HOME/libc/$ARCH/test_interp.o \
                              /$HOME/libc/$ARCH/qr_solve.o \
                              /$HOME/libc/$ARCH/r8lib.o -lm 
if [ $? -ne 0 ]; then
  echo "Errors linking and loading chebyshev_interp_1d_prb.o."
  exit
fi
#
rm chebyshev_interp_1d_prb.o
#
mv a.out chebyshev_interp_1d_prb
./chebyshev_interp_1d_prb > chebyshev_interp_1d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running chebyshev_interp_1d_prb."
  exit
fi
rm chebyshev_interp_1d_prb
#
echo "Program output written to chebyshev_interp_1d_prb_output.txt"
