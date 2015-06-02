#!/bin/bash
#
gcc -c fem_to_triangle.c
if [ $? -ne 0 ]; then
  echo "Errors compiling fem_to_triangle.c"
  exit
fi
#
gcc fem_to_triangle.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem_to_triangle.o."
  exit
fi
#
rm fem_to_triangle.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/fem_to_triangle
#
echo "Executable installed as ~/binc/$ARCH/fem_to_triangle"
