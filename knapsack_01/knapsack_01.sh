#!/bin/bash
#
cp knapsack_01.h /$HOME/include
#
gcc -c -I/$HOME/include knapsack_01.c
if [ $? -ne 0 ]; then
  echo "Errors compiling knapsack_01.c"
  exit
fi
#
mv knapsack_01.o ~/libc/$ARCH/knapsack_01.o
#
echo "Library installed as ~/libc/$ARCH/knapsack_01.o"
