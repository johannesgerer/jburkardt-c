#!/bin/bash
#
gcc -c -g -I/$HOME/include toms446_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms446_prb.c"
  exit
fi
rm compiler.txt
#
gcc toms446_prb.o /$HOME/libc/$ARCH/toms446.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms446_prb.o."
  exit
fi
#
rm toms446_prb.o
#
mv a.out toms446_prb
./toms446_prb > toms446_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms446_prb."
  exit
fi
rm toms446_prb
#
echo "Program output written to toms446_prb_output.txt"
