#!/bin/bash
#
cp subset_sum_serial.h /$HOME/include
#
gcc -c -g -I/$HOME/include subset_sum_serial.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling subset_sum_serial.c"
  exit
fi
rm compiler.txt
#
mv subset_sum_serial.o ~/libc/$ARCH/subset_sum_serial.o
#
echo "Library installed as ~/libc/$ARCH/subset_sum_serial.o"
