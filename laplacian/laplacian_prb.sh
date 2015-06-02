#!/bin/bash
#
gcc -c -g -I/$HOME/include laplacian_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling laplacian_prb.c"
  exit
fi
rm compiler.txt
#
gcc laplacian_prb.o /$HOME/libc/$ARCH/laplacian.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading laplacian_prb.o."
  exit
fi
#
rm laplacian_prb.o
#
mv a.out laplacian_prb
./laplacian_prb > laplacian_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running laplacian_prb."
  exit
fi
rm laplacian_prb
#
echo "Program output written to laplacian_prb_output.txt"
