#!/bin/bash
#
cp polygon_moments.h /$HOME/include
#
gcc -c -g -I/$HOME/include polygon_moments.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling polygon_moments.c"
  exit
fi
rm compiler.txt
#
mv polygon_moments.o ~/libc/$ARCH/polygon_moments.o
#
echo "Library installed as ~/libc/$ARCH/polygon_moments.o"
