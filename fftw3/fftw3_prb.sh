#!/bin/bash
#
gcc -c -g -I/$HOME/include fftw3_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fftw3_prb.c."
  exit
fi
rm compiler.txt
#
gcc fftw3_prb.o -lfftw3 -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fftw3_prb.o."
  exit
fi
#
rm fftw3_prb.o
#
mv a.out fftw3_prb
./fftw3_prb > fftw3_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fftw3_prb."
  exit
fi
rm fftw3_prb
#
echo "Program output written to fftw3_prb_output.txt"
