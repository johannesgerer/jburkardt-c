#!/bin/bash
#
cp circle_integrals.h /$HOME/include
#
gcc -c -g -I /$HOME/include circle_integrals.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling circle_integrals.c"
  exit
fi
rm compiler.txt
#
mv circle_integrals.o ~/libc/$ARCH/circle_integrals.o
#
echo "Library installed as ~/libc/$ARCH/circle_integrals.o"
