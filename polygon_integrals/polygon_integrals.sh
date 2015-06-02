#!/bin/bash
#
cp polygon_integrals.h /$HOME/include
#
gcc -c -g -I/$HOME/include polygon_integrals.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling polygon_integrals.c"
  exit
fi
rm compiler.txt
#
mv polygon_integrals.o ~/libc/$ARCH/polygon_integrals.o
#
echo "Library installed as ~/libc/$ARCH/polygon_integrals.o"
