#!/bin/bash
#
#  Compile and load as an executable.
#
gcc -c triangle.c
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle.c."
  exit
fi
#
gcc triangle.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangle.c."
  exit
fi
rm triangle.o
mv a.out ~/binc/$ARCH/triangle
echo "Executable installed as ~/binc/$ARCH/triangle"
#
#  Using the TRILIBRARY switch, the file can be compiled as a library.
#
gcc -c -DTRILIBRARY -g triangle.c
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle.c."
  exit
fi
#
mv triangle.o ~/libc/$ARCH/triangle.o
#
echo "Library installed as ~/libc/$ARCH/triangle.o"
