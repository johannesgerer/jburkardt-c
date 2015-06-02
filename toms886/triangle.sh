#!/bin/bash
#
gcc -c -g triangle.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle.c"
  exit
fi
rm compiler.txt
#
gcc triangle.o $HOME/libc/$ARCH/toms886.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangle.o"
  exit
fi
rm triangle.o
#
mv a.out triangle
./triangle > triangle_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running triangle"
  exit
fi
rm triangle
#
echo "Test results written to triangle_output.txt."
