#!/bin/bash
#
gcc -c -g colored_noise_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling colored_noise_prb.c."
  exit
fi
rm compiler.txt
#
gcc colored_noise_prb.o /$HOME/libc/$ARCH/colored_noise.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading colored_noise_prb.o."
  exit
fi
#
rm colored_noise_prb.o
#
mv a.out colored_noise_prb
./colored_noise_prb > colored_noise_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running colored_noise_prb."
  exit
fi
rm colored_noise_prb
#
echo "Program output written to colored_noise_prb_output.txt"
