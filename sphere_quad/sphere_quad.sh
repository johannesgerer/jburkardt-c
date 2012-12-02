#!/bin/bash
#
cp sphere_quad.h /$HOME/include
#
gcc -c -g -I /$HOME/include sphere_quad.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_quad.c."
  exit
fi
rm compiler.txt
#
mv sphere_quad.o ~/libc/$ARCH/sphere_quad.o
#
echo "Library installed as ~/libc/$ARCH/sphere_quad.o"
