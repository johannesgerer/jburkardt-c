#!/bin/bash
#
gcc -c -g unicycle_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling unicycle_prb.c"
  exit
fi
rm compiler.txt
#
gcc unicycle_prb.o /$HOME/libc/$ARCH/unicycle.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading unicycle_prb.o."
  exit
fi
#
rm unicycle_prb.o
#
mv a.out unicycle_prb
./unicycle_prb > unicycle_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running unicycle_prb."
  exit
fi
rm unicycle_prb
#
echo "Program output written to unicycle_prb_output.txt"
