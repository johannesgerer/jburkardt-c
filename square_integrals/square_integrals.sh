#!/bin/bash
#
cp square_integrals.h /$HOME/include
#
gcc -c -g -I /$HOME/include square_integrals.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling square_integrals.c"
  exit
fi
rm compiler.txt
#
mv square_integrals.o ~/libc/$ARCH/square_integrals.o
#
echo "Library installed as ~/libc/$ARCH/square_integrals.o"
