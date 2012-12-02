#!/bin/bash
#
gcc -c -g fair_dice_simulation.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fair_dice_simulation.c"
  exit
fi
rm compiler.txt
#
gcc fair_dice_simulation.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fair_dice_simulation.o"
  exit
fi
rm fair_dice_simulation.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/fair_dice_simulation
#
echo "Executable installed as ~/binc/$ARCH/fair_dice_simulation"
