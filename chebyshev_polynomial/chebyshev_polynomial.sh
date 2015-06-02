#!/bin/bash
#
cp chebyshev_polynomial.h /$HOME/include
#
gcc -c -g -I /$HOME/include chebyshev_polynomial.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling chebyshev_polynomial.c"
  exit
fi
rm compiler.txt
#
mv chebyshev_polynomial.o ~/libc/$ARCH/chebyshev_polynomial.o
#
echo "Library installed as ~/libc/$ARCH/chebyshev_polynomial.o"
