#!/bin/bash
#
gcc -c -g kronrod_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling kronrod_prb.c."
  exit
fi
rm compiler.txt
#
gcc kronrod_prb.o /$HOME/libc/$ARCH/kronrod.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading kronrod_prb.o."
  exit
fi
#
rm kronrod_prb.o
#
mv a.out kronrod_prb
./kronrod_prb > kronrod_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running kronrod_prb."
  exit
fi
rm kronrod_prb
#
echo "Program output written to kronrod_prb_output.txt"
