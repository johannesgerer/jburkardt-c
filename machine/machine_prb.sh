#!/bin/bash
#
gcc -c -g machine_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling machine_prb.c."
  exit
fi
rm compiler.txt
#
gcc machine_prb.o /$HOME/libc/$ARCH/machine.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading machine_prb.o."
  exit
fi
#
rm machine_prb.o
#
mv a.out machine_prb
./machine_prb > machine_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running machine_prb."
  exit
fi
rm machine_prb
#
echo "Program output written to machine_prb_output.txt"
