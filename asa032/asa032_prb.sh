#!/bin/bash
#
gcc -c -g asa032_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa032_prb.c."
  exit
fi
rm compiler.txt
#
gcc asa032_prb.o /$HOME/libc/$ARCH/asa032.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa032_prb.o."
  exit
fi
#
rm asa032_prb.o
#
mv a.out asa032_prb
./asa032_prb > asa032_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa032_prb."
  exit
fi
rm asa032_prb
#
echo "Program output written to asa032_prb_output.txt"
