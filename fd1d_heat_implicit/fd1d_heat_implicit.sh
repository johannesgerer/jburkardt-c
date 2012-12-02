#!/bin/bash
#
gcc -c fd1d_heat_implicit.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fd1d_heat_implicit.c"
  exit
fi
rm compiler.txt
#
gcc fd1d_heat_implicit.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fd1d_heat_implicit.o"
  exit
fi
rm fd1d_heat_implicit.o
#
mv a.out ~/binc/$ARCH/fd1d_heat_implicit
#
echo "Executable installed as ~/binc/$ARCH/fd1d_heat_implicit"
