#!/bin/bash
#
gcc -c -g -I/$HOME/include haar_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling haar_prb.c."
  exit
fi
rm compiler.txt
#
gcc haar_prb.o /$HOME/libc/$ARCH/haar.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading haar_prb.o."
  exit
fi
#
rm haar_prb.o
#
mv a.out haar_prb
./haar_prb > haar_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running haar_prb."
  exit
fi
rm haar_prb
#
echo "Program output written to haar_prb_output.txt"
