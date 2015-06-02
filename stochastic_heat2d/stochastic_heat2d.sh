#!/bin/bash
#
cp stochastic_heat2d.h /$HOME/include
#
gcc -c -g -I/$HOME/include stochastic_heat2d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling stochastic_heat2d.c"
  exit
fi
rm compiler.txt
#
mv stochastic_heat2d.o ~/libc/$ARCH/stochastic_heat2d.o
#
echo "Library installed as ~/libc/$ARCH/stochastic_heat2d.o"
