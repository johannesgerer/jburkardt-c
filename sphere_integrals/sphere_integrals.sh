#!/bin/bash
#
cp sphere_integrals.h /$HOME/include
#
gcc -c -g -I/$HOME/include sphere_integrals.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_integrals.c"
  exit
fi
rm compiler.txt
#
mv sphere_integrals.o ~/libc/$ARCH/sphere_integrals.o
#
echo "Library installed as ~/libc/$ARCH/sphere_integrals.o"
