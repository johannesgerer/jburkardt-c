#!/bin/bash
#
gcc -c -g fem2d_poisson_rectangle.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_poisson_rectangle.c"
  exit
fi
rm compiler.txt
#
gcc fem2d_poisson_rectangle.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem2d_poisson_rectangle.o"
  exit
fi
#
rm fem2d_poisson_rectangle.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/fem2d_poisson_rectangle
#
echo "Executable installed as ~/binc/$ARCH/fem2d_poisson_rectangle"
