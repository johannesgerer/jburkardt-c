#!/bin/bash
#
gcc -c -g -I/$HOME/include chebyshev_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling chebyshev_prb.c"
  exit
fi
rm compiler.txt
#
gcc chebyshev_prb.o /$HOME/libc/$ARCH/chebyshev.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading chebyshev_prb.o."
  exit
fi
#
rm chebyshev_prb.o
#
mv a.out chebyshev_prb
./chebyshev_prb > chebyshev_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running chebyshev_prb."
  exit
fi
rm chebyshev_prb
#
echo "Program output written to chebyshev_prb_output.txt"
