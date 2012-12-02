#!/bin/bash
#
gcc -c forest_fire_simulation.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling forest_fire_simulation.c"
  exit
fi
rm compiler.txt
#
gcc -O -L/usr/X11R6/lib forest_fire_simulation.c -lX11 
if [ $? -ne 0 ]; then
  echo "Errors linking and loading forest_fire_simulation.o"
  exit
fi
#
rm forest_fire_simulation.o
mv a.out ~/binc/$ARCH/forest_fire_simulation
#
echo "Program installed as ~/binc/$ARCH/forest_fire_simulation"
