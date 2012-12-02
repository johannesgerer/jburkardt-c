#!/bin/bash
#
gcc -c string_simulation.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling string_simulation.c"
  exit
fi
rm compiler.txt
#
gcc string_simulation.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking string_simulation.o"
  exit
fi
#
rm string_simulation.o
mv a.out ~/binc/$ARCH/string_simulation
#
echo "Executable installed as ~/binc/$ARCH/string_simulation"
