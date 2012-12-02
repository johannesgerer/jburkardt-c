#!/bin/bash
#
gcc -c -g asa239_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa239_prb.c."
  exit
fi
rm compiler.txt
#
gcc asa239_prb.o /$HOME/libc/$ARCH/asa239.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa239_prb.o."
  exit
fi
#
rm asa239_prb.o
#
mv a.out asa239_prb
./asa239_prb > asa239_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa239_prb."
  exit
fi
rm asa239_prb
#
echo "Program output written to asa239_prb_output.txt"
