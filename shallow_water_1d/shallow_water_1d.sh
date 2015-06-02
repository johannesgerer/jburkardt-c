#!/bin/bash
#
gcc -c -g shallow_water_1d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling shallow_water_1d.c"
  exit
fi
rm compiler.txt
#
gcc shallow_water_1d.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading shallow_water_1d.o"
  exit
fi
rm shallow_water_1d.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/shallow_water_1d
#
echo "Program installed as ~/binc/$ARCH/shallow_water_1d"
