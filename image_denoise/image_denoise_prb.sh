#!/bin/bash
#
gcc -c -g -I/$HOME/include image_denoise_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling image_denoise_prb.c"
  exit
fi
rm compiler.txt
#
gcc image_denoise_prb.o /$HOME/libc/$ARCH/image_denoise.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading image_denoise_prb.o."
  exit
fi
#
rm image_denoise_prb.o
#
mv a.out image_denoise_prb
./image_denoise_prb > image_denoise_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running image_denoise_prb."
  exit
fi
rm image_denoise_prb
#
echo "Program output written to image_denoise_prb_output.txt"
