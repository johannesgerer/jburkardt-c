#!/bin/bash
#
cp chebyshev.h /$HOME/include
#
gcc -c -g -I /$HOME/include chebyshev.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling chebyshev.c"
  exit
fi
rm compiler.txt
#
mv chebyshev.o ~/libc/$ARCH/chebyshev.o
#
echo "Library installed as ~/libc/$ARCH/chebyshev.o"
