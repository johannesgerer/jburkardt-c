#!/bin/bash
#
cp subset_sum.h /$HOME/include
#
gcc -c -g -I /$HOME/include subset_sum.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling subset_sum.c"
  exit
fi
rm compiler.txt
#
mv subset_sum.o ~/libc/$ARCH/subset_sum.o
#
echo "Library installed as ~/libc/$ARCH/subset_sum.o"
