#!/bin/bash
#
cp wavelet.h /$HOME/include
#
gcc -c -g -I /$HOME/include wavelet.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling wavelet.c"
  exit
fi
rm compiler.txt
#
mv wavelet.o ~/libc/$ARCH/wavelet.o
#
echo "Library installed as ~/libc/$ARCH/wavelet.o"
