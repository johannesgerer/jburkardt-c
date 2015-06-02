#!/bin/bash
#
cp stochastic_diffusion.h /$HOME/include
#
gcc -c -g -I/$HOME/include stochastic_diffusion.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling stochastic_diffusion.c"
  exit
fi
rm compiler.txt
#
mv stochastic_diffusion.o ~/libc/$ARCH/stochastic_diffusion.o
#
echo "Library installed as ~/libc/$ARCH/stochastic_diffusion.o"
