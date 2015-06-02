#!/bin/bash
#
cp chebyshev_series.h /$HOME/include
#
gcc -c -I /$HOME/include chebyshev_series.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling chebyshev_series.c"
  exit
fi
rm compiler.txt
#
mv chebyshev_series.o ~/libc/$ARCH/chebyshev_series.o
#
echo "Library installed as ~/libc/$ARCH/chebyshev_series.o"
