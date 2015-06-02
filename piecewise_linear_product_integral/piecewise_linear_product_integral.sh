#!/bin/bash
#
cp piecewise_linear_product_integral.h /$HOME/include
#
gcc -c -g -I /$HOME/include piecewise_linear_product_integral.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling piecewise_linear_product_integral.c"
  exit
fi
rm compiler.txt
#
mv piecewise_linear_product_integral.o ~/libc/$ARCH/piecewise_linear_product_integral.o
#
echo "Library installed as ~/libc/$ARCH/piecewise_linear_product_integral.o"
