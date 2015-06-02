#!/bin/bash
#
gcc -c wedge_exactness.c
if [ $? -ne 0 ]; then
  echo "Errors compiling wedge_exactness.c"
  exit
fi
#
gcc wedge_exactness.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading wedge_exactness.o"
  exit
fi
rm wedge_exactness.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/wedge_exactness
#
echo "Executable installed as ~/binc/$ARCH/wedge_exactness"
