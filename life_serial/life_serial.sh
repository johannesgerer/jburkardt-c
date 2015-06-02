#!/bin/bash
#
gcc -c -g life_serial.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling life_serial.c."
  exit
fi
rm compiler.txt
#
gcc life_serial.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading life_serial.o."
  exit
fi
#
rm life_serial.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/life_serial
#
echo "Executable installed as ~/binc/$ARCH/life_serial"
