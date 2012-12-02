#!/bin/bash
#
gcc -c -g brent_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling brent_prb.c."
  exit
fi
rm compiler.txt
#
gcc brent_prb.o /$HOME/libc/$ARCH/brent.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading brent_prb.o."
  exit
fi
#
rm brent_prb.o
#
mv a.out brent_prb
./brent_prb > brent_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running brent_prb."
  exit
fi
rm brent_prb
#
echo "Program output written to brent_prb_output.txt"
