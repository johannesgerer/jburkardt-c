#!/bin/bash
#
gcc -c -g pbmlib_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pbmlib_prb.c."
  exit
fi
rm compiler.txt
#
gcc pbmlib_prb.o /$HOME/libc/$ARCH/pbmlib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pbmlib_prb.o."
  exit
fi
#
rm pbmlib_prb.o
#
mv a.out pbmlib_prb
./pbmlib_prb > pbmlib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pbmlib_prb."
  exit
fi
rm pbmlib_prb
#
echo "Program output written to pbmlib_prb_output.txt"
