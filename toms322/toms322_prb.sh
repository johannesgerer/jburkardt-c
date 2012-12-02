#!/bin/bash
#
gcc -c -g toms322_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms322_prb.c."
  exit
fi
rm compiler.txt
#
gcc toms322_prb.o /$HOME/libc/$ARCH/toms322.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms322_prb.o."
  exit
fi
#
rm toms322_prb.o
#
mv a.out toms322_prb
./toms322_prb > toms322_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms322_prb."
  exit
fi
rm toms322_prb
#
echo "Program output written to toms322_prb_output.txt"
