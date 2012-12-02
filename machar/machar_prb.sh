#!/bin/bash
#
gcc -c -g machar_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling machar_prb.c."
  exit
fi
rm compiler.txt
#
gcc machar_prb.o /$HOME/libc/$ARCH/machar.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading machar_prb.o."
  exit
fi
#
rm machar_prb.o
#
mv a.out machar_prb
./machar_prb > machar_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running machar_prb."
  exit
fi
rm machar_prb
#
echo "Program output written to machar_prb_output.txt"
