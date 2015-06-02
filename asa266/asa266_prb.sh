#!/bin/bash
#
gcc -c -g asa266_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa266_prb.c."
  exit
fi
rm compiler.txt
#
gcc asa266_prb.o /$HOME/libc/$ARCH/asa266.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa266_prb.o."
  exit
fi
#
rm asa266_prb.o
#
mv a.out asa266_prb
./asa266_prb > asa266_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa266_prb."
  exit
fi
rm asa266_prb
#
echo "Program output written to asa266_prb_output.txt"
