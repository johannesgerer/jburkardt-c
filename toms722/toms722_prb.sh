#!/bin/bash
#
gcc -c -g toms722_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms722_prb.c."
  exit
fi
rm compiler.txt
#
gcc toms722_prb.o /$HOME/libc/$ARCH/toms722.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms722_prb.o."
  exit
fi
#
rm toms722_prb.o
#
mv a.out toms722_prb
./toms722_prb > toms722_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms722_prb."
  exit
fi
rm toms722_prb
#
echo "Program output written to toms722_prb_output.txt"
