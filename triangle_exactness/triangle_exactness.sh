#!/bin/bash
#
gcc -c -I$HOME/include triangle_exactness.c
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_exactness.c"
  exit
fi
#
gcc triangle_exactness.o -lm 
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangle_exactness.o."
  exit
fi
#
rm triangle_exactness.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/triangle_exactness
#
echo "Executable installed as ~/binc/$ARCH/triangle_exactness"
