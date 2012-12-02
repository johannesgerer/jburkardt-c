#!/bin/bash
#
gcc -c -g pink_noise_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pink_noise_prb.c."
  exit
fi
rm compiler.txt
#
gcc pink_noise_prb.o /$HOME/libc/$ARCH/pink_noise.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pink_noise_prb.o."
  exit
fi
#
rm pink_noise_prb.o
#
mv a.out pink_noise_prb
./pink_noise_prb > pink_noise_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pink_noise_prb."
  exit
fi
rm pink_noise_prb
#
echo "Program output written to pink_noise_prb_output.txt"
