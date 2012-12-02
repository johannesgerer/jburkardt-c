#!/bin/bash
#
gcc -c -g asa005_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa005_prb.c."
  exit
fi
rm compiler.txt
#
gcc asa005_prb.o /$HOME/libc/$ARCH/asa005.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa005_prb.o."
  exit
fi
#
rm asa005_prb.o
#
mv a.out asa005_prb
./asa005_prb > asa005_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa005_prb."
  exit
fi
rm asa005_prb
#
echo "Program output written to asa005_prb_output.txt"
