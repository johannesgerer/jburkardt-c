#!/bin/bash
#
gcc -c triangle_to_medit.c
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_to_medit.c"
  exit
fi
#
gcc triangle_to_medit.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangle_to_medit.o."
  exit
fi
rm triangle_to_medit.o
#
mv a.out ~/binc/$ARCH/triangle_to_medit
#
echo "Executable installed as ~/binc/$ARCH/triangle_to_medit"
