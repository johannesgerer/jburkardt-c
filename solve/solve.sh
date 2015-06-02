#!/bin/bash
#
cp solve.h /$HOME/include
#
gcc -c -I/$HOME/include solve.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling solve.c"
  exit
fi
rm compiler.txt
#
mv solve.o ~/libc/$ARCH/solve.o
#
echo "Library installed as ~/libc/$ARCH/solve.o"
