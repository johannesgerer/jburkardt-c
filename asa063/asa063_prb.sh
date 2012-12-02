#!/bin/bash
#
gcc -c -g asa063_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa063_prb.c."
  exit
fi
rm compiler.txt
#
gcc asa063_prb.o /$HOME/libc/$ARCH/asa063.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa063_prb.o."
  exit
fi
#
rm asa063_prb.o
#
mv a.out asa063_prb
./asa063_prb > asa063_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa063_prb."
  exit
fi
rm asa063_prb
#
echo "Program output written to asa063_prb_output.txt"
