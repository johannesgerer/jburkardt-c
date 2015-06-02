#!/bin/bash
#
gcc -c -g -I/$HOME/include ou_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ou_prb.c"
  exit
fi
rm compiler.txt
#
gcc ou_prb.o /$HOME/libc/$ARCH/ou.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ou_prb.o."
  exit
fi
#
rm ou_prb.o
#
mv a.out ou_prb
./ou_prb > ou_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ou_prb."
  exit
fi
rm ou_prb
#
echo "Program output written to ou_prb_output.txt"
