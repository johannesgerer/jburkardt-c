#!/bin/bash
#
gcc -c -g -I/$HOME/include cycle_brent_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cycle_brent_prb.c"
  exit
fi
rm compiler.txt
#
gcc cycle_brent_prb.o /$HOME/libc/$ARCH/cycle_brent.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cycle_brent_prb.o"
  exit
fi
#
rm cycle_brent_prb.o
#
mv a.out cycle_brent_prb
./cycle_brent_prb > cycle_brent_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cycle_brent_prb."
  exit
fi
rm cycle_brent_prb
#
echo "Program output written to cycle_brent_prb_output.txt"
