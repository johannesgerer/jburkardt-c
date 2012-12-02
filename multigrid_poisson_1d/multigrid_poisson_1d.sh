#!/bin/bash
#
cp multigrid_poisson_1d.h /$HOME/include
#
gcc -c -g multigrid_poisson_1d.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling multigrid_poisson_1d.c."
  exit
fi
rm compiler.txt
#
mv multigrid_poisson_1d.o ~/libc/$ARCH/multigrid_poisson_1d.o
#
echo "Library installed as ~/libc/$ARCH/multigrid_poisson_1d.o"
