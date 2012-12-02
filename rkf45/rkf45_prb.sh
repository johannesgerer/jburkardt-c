#!/bin/bash
#
gcc -c -g rkf45_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rkf45_prb.c."
  exit
fi
rm compiler.txt
#
gcc rkf45_prb.o /$HOME/libc/$ARCH/rkf45.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading rkf45_prb.o."
  exit
fi
#
rm rkf45_prb.o
#
mv a.out rkf45_prb
./rkf45_prb > rkf45_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running rkf45_prb."
  exit
fi
rm rkf45_prb
#
echo "Program output written to rkf45_prb_output.txt"
