#!/bin/bash
#
gcc -c -g fem2d_poisson_rectangle_linear.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_poisson_rectangle_linear.c"
  exit
fi
rm compiler.txt
#
gcc fem2d_poisson_rectangle_linear.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem2d_poisson_rectangle_linear.o"
  exit
fi
#
rm fem2d_poisson_rectangle_linear.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/fem2d_poisson_rectangle_linear
#
echo "Executable installed as ~/binc/$ARCH/fem2d_poisson_rectangle_linear"
