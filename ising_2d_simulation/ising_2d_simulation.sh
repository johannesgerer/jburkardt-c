#!/bin/bash
#
gcc -c -g ising_2d_simulation.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ising_2d_simulation.c."
  exit
fi
rm compiler.txt
#
gcc ising_2d_simulation.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ising_2d_simulation.o."
  exit
fi
#
rm ising_2d_simulation.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/ising_2d_simulation
#
echo "Executable installed as ~/binc/$ARCH/ising_2d_simulation"
