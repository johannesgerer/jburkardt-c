#!/bin/bash
#
gcc -c fd1d_advection_lax_wendroff.c
if [ $? -ne 0 ]; then
  echo "Errors compiling fd1d_advection_lax_wendroff.c"
  exit
fi
#
gcc fd1d_advection_lax_wendroff.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking fd1d_advection_lax_wendroff.o"
  exit
fi
#
rm fd1d_advection_lax_wendroff.o
mv a.out ~/binc/$ARCH/fd1d_advection_lax_wendroff
#
echo "Executable installed as ~/binc/$ARCH/fd1d_advection_lax_wendroff"
