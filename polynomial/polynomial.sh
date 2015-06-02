#!/bin/bash
#
cp polynomial.h /$HOME/include
#
gcc -c -I/$HOME/include polynomial.c
if [ $? -ne 0 ]; then
  echo "Errors compiling polynomial.c"
  exit
fi
#
mv polynomial.o ~/libc/$ARCH/polynomial.o
#
echo "Library installed as ~/libc/$ARCH/polynomial.o"
