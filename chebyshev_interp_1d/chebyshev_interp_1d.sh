#!/bin/bash
#
cp chebyshev_interp_1d.h /$HOME/include
#
gcc -c -g -I/$HOME/include chebyshev_interp_1d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling chebyshev_interp_1d.c"
  exit
fi
rm compiler.txt
#
mv chebyshev_interp_1d.o ~/libc/$ARCH/chebyshev_interp_1d.o
#
echo "Library installed as ~/libc/$ARCH/chebyshev_interp_1d.o"
