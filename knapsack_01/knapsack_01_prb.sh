#!/bin/bash
#
gcc -c -I/$HOME/include knapsack_01_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling knapsack_01_prb.c"
  exit
fi
#
gcc -o knapsack_01_prb knapsack_01_prb.o /$HOME/libc/$ARCH/knapsack_01.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading knapsack_01_prb.o."
  exit
fi
#
rm knapsack_01_prb.o
#
./knapsack_01_prb > knapsack_01_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running knapsack_01_prb."
  exit
fi
rm knapsack_01_prb
#
echo "Program output written to knapsack_01_prb_output.txt"
