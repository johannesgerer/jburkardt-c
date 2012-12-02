#!/bin/bash
#
gcc -c -g asa159_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa159_prb.c."
  exit
fi
rm compiler.txt
#
gcc asa159_prb.o /$HOME/libc/$ARCH/asa159.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa159_prb.o."
  exit
fi
#
rm asa159_prb.o
#
mv a.out asa159_prb
./asa159_prb > asa159_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa159_prb."
  exit
fi
rm asa159_prb
#
echo "Program output written to asa159_prb_output.txt"
