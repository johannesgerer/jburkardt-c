#!/bin/bash
#
cp triangle_svg.h /$HOME/include
#
gcc -c -I/$HOME/include triangle_svg.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_svg.c"
  exit
fi
rm compiler.txt
#
mv triangle_svg.o ~/libc/$ARCH/triangle_svg.o
#
echo "Library installed as ~/libc/$ARCH/triangle_svg.o"
