#!/bin/bash
#
gcc -c -I$HOME/include hypercube_exactness.c
if [ $? -ne 0 ]; then
  echo "Errors compiling hypercube_exactness.c"
  exit
fi
#
gcc hypercube_exactness.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hypercube_exactness.o."
  exit
fi
#
rm hypercube_exactness.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/hypercube_exactness
#
echo "Executable installed as ~/binc/$ARCH/hypercube_exactness"
