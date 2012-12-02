#!/bin/bash
#
gcc -c -g asa121_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa121_prb.c."
  exit
fi
rm compiler.txt
#
gcc asa121_prb.o /$HOME/libc/$ARCH/asa121.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa121_prb.o."
  exit
fi
#
rm asa121_prb.o
#
mv a.out asa121_prb
./asa121_prb > asa121_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa121_prb."
  exit
fi
rm asa121_prb
#
echo "Program output written to asa121_prb_output.txt"
