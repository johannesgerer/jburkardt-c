#!/bin/bash
#
gcc -c fd1d_advection_diffusion_steady.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fd1d_advection_diffusion_steady.c"
  exit
fi
rm compiler.txt
#
gcc fd1d_advection_diffusion_steady.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking fd1d_advection_diffusion_steady.o"
  exit
fi
#
rm fd1d_advection_diffusion_steady.o
mv a.out ~/binc/$ARCH/fd1d_advection_diffusion_steady
#
echo "Executable installed as ~/binc/$ARCH/fd1d_advection_diffusion_steady"
