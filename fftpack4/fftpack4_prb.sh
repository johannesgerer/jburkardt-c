#!/bin/bash
#
gcc -c -g -I/$HOME/include fftpack4_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fftpack4_prb.c"
  exit
fi
rm compiler.txt
#
gcc fftpack4_prb.o /$HOME/libc/$ARCH/fftpack4.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fftpack4_prb.o"
  exit
fi
#
rm fftpack4_prb.o
#
mv a.out fftpack4_prb
./fftpack4_prb > fftpack4_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fftpack4_prb."
  exit
fi
rm fftpack4_prb
#
echo "Program output written to fftpack4_prb_output.txt"
