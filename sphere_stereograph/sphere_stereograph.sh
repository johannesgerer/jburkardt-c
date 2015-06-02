#!/bin/bash
#
cp sphere_stereograph.h /$HOME/include
#
gcc -c -g -I /$HOME/include sphere_stereograph.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_stereograph.c"
  exit
fi
rm compiler.txt
#
mv sphere_stereograph.o ~/libc/$ARCH/sphere_stereograph.o
#
echo "Library installed as ~/libc/$ARCH/sphere_stereograph.o"
