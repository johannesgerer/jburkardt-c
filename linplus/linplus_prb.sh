#!/bin/bash
#
gcc -c linplus_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling linplus_prb.c."
  exit
fi
#
gcc linplus_prb.o /$HOME/libc/$ARCH/linplus.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading linplus_prb.o."
  exit
fi
#
rm linplus_prb.o
#
mv a.out linplus_prb
./linplus_prb > linplus_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running linplus_prb."
  exit
fi
rm linplus_prb
#
echo "Program output written to linplus_prb_output.txt"
