#!/bin/bash
#
gcc -c medit_to_fem.c
if [ $? -ne 0 ]; then
  echo "Errors compiling medit_to_fem.c"
  exit
fi
#
gcc medit_to_fem.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading medit_to_fem.o."
  exit
fi
#
rm medit_to_fem.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/medit_to_fem
#
echo "Executable installed as ~/binc/$ARCH/medit_to_fem"
