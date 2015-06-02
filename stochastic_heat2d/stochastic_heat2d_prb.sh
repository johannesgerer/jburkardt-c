#!/bin/bash
#
gcc -c -g -I/$HOME/include stochastic_heat2d_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling stochastic_heat2d_prb.c"
  exit
fi
rm compiler.txt
#
gcc stochastic_heat2d_prb.o /$HOME/libc/$ARCH/stochastic_heat2d.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading stochastic_heat2d_prb.o."
  exit
fi
#
rm stochastic_heat2d_prb.o
#
mv a.out stochastic_heat2d_prb
./stochastic_heat2d_prb > stochastic_heat2d_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running stochastic_heat2d_prb."
  exit
fi
rm stochastic_heat2d_prb
#
echo "Program output written to stochastic_heat2d_prb_output.txt"
