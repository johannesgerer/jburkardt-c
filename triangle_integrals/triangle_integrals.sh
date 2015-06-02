#!/bin/bash
#
cp triangle_integrals.h /$HOME/include
#
gcc -c -I/$HOME/include triangle_integrals.c
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_integrals.c"
  exit
fi
#
mv triangle_integrals.o ~/libc/$ARCH/triangle_integrals.o
#
echo "Library installed as ~/libc/$ARCH/triangle_integrals.o"
