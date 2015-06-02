#!/bin/bash
#
cp simplex_integrals.h /$HOME/include
#
gcc -c -g -I /$HOME/include simplex_integrals.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling simplex_integrals.c"
  exit
fi
rm compiler.txt
#
mv simplex_integrals.o ~/libc/$ARCH/simplex_integrals.o
#
echo "Library installed as ~/libc/$ARCH/simplex_integrals.o"
