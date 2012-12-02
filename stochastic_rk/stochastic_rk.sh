#!/bin/bash
#
cp stochastic_rk.h /$HOME/include
#
gcc -c -g stochastic_rk.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling stochastic_rk.c."
  exit
fi
rm compiler.txt
#
mv stochastic_rk.o ~/libc/$ARCH/stochastic_rk.o
#
echo "Library installed as ~/libc/$ARCH/stochastic_rk.o"
