#!/bin/bash
#
gcc -c -g asa058_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling asa058_prb.c."
  exit
fi
rm compiler.txt
#
gcc asa058_prb.o /$HOME/libc/$ARCH/asa058.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading asa058_prb.o."
  exit
fi
#
rm asa058_prb.o
#
mv a.out asa058_prb
./asa058_prb > asa058_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running asa058_prb."
  exit
fi
rm asa058_prb
#
echo "Program output written to asa058_prb_output.txt"
