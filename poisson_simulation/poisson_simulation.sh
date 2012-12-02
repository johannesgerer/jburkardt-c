#!/bin/bash
#
cp poisson_simulation.h /$HOME/include
#
gcc -c -g -I /$HOME/include poisson_simulation.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling poisson_simulation.c"
  exit
fi
rm compiler.txt
#
mv poisson_simulation.o ~/libc/$ARCH/poisson_simulation.o
#
echo "Library installed as ~/libc/$ARCH/poisson_simulation.o"
