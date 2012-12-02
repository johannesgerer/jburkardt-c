#!/bin/bash
#
cp triangulation.h /$HOME/include
#
gcc -c -g -I /$HOME/include triangulation.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation.c"
  exit
fi
rm compiler.txt
#
mv triangulation.o ~/libc/$ARCH/triangulation.o
#
echo "Library installed as ~/libc/$ARCH/triangulation.o"
