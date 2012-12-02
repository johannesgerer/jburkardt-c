#!/bin/bash
#
gcc -c -g owens_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling owens_prb.c."
  exit
fi
rm compiler.txt
#
gcc owens_prb.o /$HOME/libc/$ARCH/owens.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading owens_prb.o."
  exit
fi
#
rm owens_prb.o
#
mv a.out owens_prb
./owens_prb > owens_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running owens_prb."
  exit
fi
rm owens_prb
#
echo "Program output written to owens_prb_output.txt"
