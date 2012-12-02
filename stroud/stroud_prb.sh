#!/bin/bash
#
gcc -c -g stroud_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling stroud_prb.c."
  exit
fi
rm compiler.txt
#
gcc stroud_prb.o /$HOME/libc/$ARCH/stroud.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading stroud_prb.o."
  exit
fi
#
rm stroud_prb.o
#
mv a.out stroud_prb
./stroud_prb > stroud_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running stroud_prb."
  exit
fi
rm stroud_prb
#
echo "Program output written to stroud_prb_output.txt"
