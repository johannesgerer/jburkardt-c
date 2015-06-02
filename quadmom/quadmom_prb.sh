#!/bin/bash
#
gcc -c -g -I/$HOME/include quadmom_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quadmom_prb.c"
  exit
fi
rm compiler.txt
#
gcc quadmom_prb.o /$HOME/libc/$ARCH/quadmom.o  /$HOME/libc/$ARCH/toms655.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quadmom_prb.o."
  exit
fi
#
rm quadmom_prb.o
#
mv a.out quadmom_prb
./quadmom_prb > quadmom_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running quadmom_prb."
  exit
fi
rm quadmom_prb
#
echo "Program output written to quadmom_prb_output.txt"
