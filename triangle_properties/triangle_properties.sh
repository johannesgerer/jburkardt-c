#!/bin/bash
#
gcc -c -g triangle_properties.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_properties.c"
  exit
fi
rm compiler.txt
#
gcc triangle_properties.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangle_properties.o"
  exit
fi
rm triangle_properties.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/triangle_properties
#
echo "Executable installed as ~/binc/$ARCH/triangle_properties"
