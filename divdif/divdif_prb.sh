#!/bin/bash
#
gcc -c -g -I/$HOME/include divdif_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling divdif_prb.c."
  exit
fi
rm compiler.txt
#
gcc divdif_prb.o /$HOME/libc/$ARCH/divdif.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading divdif_prb.o."
  exit
fi
#
rm divdif_prb.o
#
mv a.out divdif_prb
./divdif_prb > divdif_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running divdif_prb."
  exit
fi
rm divdif_prb
#
echo "Program output written to divdif_prb_output.txt"
