#!/bin/bash
#
gcc -c -g toms179_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms179_prb.c."
  exit
fi
rm compiler.txt
#
gcc toms179_prb.o /$HOME/libc/$ARCH/toms179.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms179_prb.o."
  exit
fi
#
rm toms179_prb.o
#
mv a.out toms179_prb
./toms179_prb > toms179_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms179_prb."
  exit
fi
rm toms179_prb
#
echo "Program output written to toms179_prb_output.txt"
