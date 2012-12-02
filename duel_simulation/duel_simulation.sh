#!/bin/bash
#
gcc -c -g duel_simulation.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling duel_simulation.c"
  exit
fi
rm compiler.txt
#
gcc duel_simulation.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading duel_simulation.o"
  exit
fi
rm duel_simulation.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/duel_simulation
#
echo "Executable installed as ~/binc/$ARCH/duel_simulation"
