#!/bin/bash
#
gcc -c -g asa299_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa299_prb.c."
  exit
fi
rm compiler.txt
#
gcc asa299_prb.o /$HOME/libc/$ARCH/asa299.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa299_prb.o."
  exit
fi
#
rm asa299_prb.o
#
mv a.out asa299_prb
./asa299_prb > asa299_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa299_prb."
  exit
fi
rm asa299_prb
#
echo "Program output written to asa299_prb_output.txt"
