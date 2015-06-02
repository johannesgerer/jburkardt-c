#!/bin/bash
#
gcc -c -g -I/$HOME/include disk_integrals_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling disk_integrals_prb.c"
  exit
fi
rm compiler.txt
#
gcc disk_integrals_prb.o /$HOME/libc/$ARCH/disk_integrals.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading disk_integrals_prb.o."
  exit
fi
#
rm disk_integrals_prb.o
#
mv a.out disk_integrals_prb
./disk_integrals_prb > disk_integrals_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running disk_integrals_prb."
  exit
fi
rm disk_integrals_prb
#
echo "Program output written to disk_integrals_prb_output.txt"
