#!/bin/bash
#
cp multigrid_poisson_1d.h /$HOME/include
#
gcc -c multigrid_poisson_1d.c
if [ $? -ne 0 ]; then
  echo "Errors compiling multigrid_poisson_1d.c."
  exit
fi
#
mv multigrid_poisson_1d.o ~/libc/$ARCH/multigrid_poisson_1d.o
#
echo "Library installed as ~/libc/$ARCH/multigrid_poisson_1d.o"
