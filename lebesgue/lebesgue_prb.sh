#!/bin/bash
#
gcc -c -g -I/$HOME/include lebesgue_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lebesgue_prb.c"
  exit
fi
rm compiler.txt
#
gcc lebesgue_prb.o /$HOME/libc/$ARCH/lebesgue.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lebesgue_prb.o."
  exit
fi
#
rm lebesgue_prb.o
#
mv a.out lebesgue_prb
./lebesgue_prb > lebesgue_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lebesgue_prb."
  exit
fi
rm lebesgue_prb
#
echo "Program output written to lebesgue_prb_output.txt"
