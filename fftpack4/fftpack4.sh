#!/bin/bash
#
cp fftpack4.h /$HOME/include
cp fftpack4_precision.h /$HOME/include
#
gcc -c -g -I /$HOME/include fftpack4.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fftpack4.c"
  exit
fi
rm compiler.txt
#
mv fftpack4.o ~/libc/$ARCH/fftpack4.o
#
echo "Library installed as ~/libc/$ARCH/fftpack4.o"
