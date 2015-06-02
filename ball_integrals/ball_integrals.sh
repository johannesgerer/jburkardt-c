#!/bin/bash
#
cp ball_integrals.h /$HOME/include
#
gcc -c -g -I/$HOME/include ball_integrals.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ball_integrals.c"
  exit
fi
rm compiler.txt
#
mv ball_integrals.o ~/libc/$ARCH/ball_integrals.o
#
echo "Library installed as ~/libc/$ARCH/ball_integrals.o"
