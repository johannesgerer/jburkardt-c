#!/bin/bash
#
gcc -c -I/$HOME/include chebyshev_series_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling chebyshev_series_prb.c"
  exit
fi
rm compiler.txt
#
gcc chebyshev_series_prb.o /$HOME/libc/$ARCH/chebyshev_series.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading chebyshev_series_prb.o"
  exit
fi
#
rm chebyshev_series_prb.o
#
mv a.out chebyshev_series_prb
./chebyshev_series_prb > chebyshev_series_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running chebyshev_series_prb."
  exit
fi
rm chebyshev_series_prb
#
echo "Program output written to chebyshev_series_prb_output.txt"
