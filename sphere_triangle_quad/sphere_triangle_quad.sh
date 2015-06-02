#!/bin/bash
#
cp sphere_triangle_quad.h /$HOME/include
#
gcc -c -I /$HOME/include sphere_triangle_quad.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_triangle_quad.c"
  exit
fi
rm compiler.txt
#
mv sphere_triangle_quad.o ~/libc/$ARCH/sphere_triangle_quad.o
#
echo "Library installed as ~/libc/$ARCH/sphere_triangle_quad.o"
