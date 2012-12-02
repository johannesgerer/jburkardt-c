#!/bin/bash
#
gcc -c heated_plate.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling heated_plate.c"
  exit
fi
rm compiler.txt
#
gcc heated_plate.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking heated_plate.o"
  exit
fi
#
rm heated_plate.o
mv a.out ~/binc/$ARCH/heated_plate
#
echo "Executable installed as ~/binc/$ARCH/heated_plate"
