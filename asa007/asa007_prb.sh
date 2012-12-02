#!/bin/bash
#
gcc -c -g asa007_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa007_prb.c."
  exit
fi
rm compiler.txt
#
gcc asa007_prb.o /$HOME/libc/$ARCH/asa007.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa007_prb.o."
  exit
fi
#
rm asa007_prb.o
#
mv a.out asa007_prb
./asa007_prb > asa007_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa007_prb."
  exit
fi
rm asa007_prb
#
echo "Program output written to asa007_prb_output.txt"
