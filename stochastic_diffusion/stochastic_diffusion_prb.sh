#!/bin/bash
#
gcc -c -g -I/$HOME/include stochastic_diffusion_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling stochastic_diffusion_prb.c"
  exit
fi
rm compiler.txt
#
gcc stochastic_diffusion_prb.o /$HOME/libc/$ARCH/stochastic_diffusion.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading stochastic_diffusion_prb.o."
  exit
fi
#
rm stochastic_diffusion_prb.o
#
mv a.out stochastic_diffusion_prb
./stochastic_diffusion_prb > stochastic_diffusion_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running stochastic_diffusion_prb."
  exit
fi
rm stochastic_diffusion_prb
#
echo "Program output written to stochastic_diffusion_prb_output.txt"
