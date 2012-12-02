#!/bin/bash
#
cp collatz_recursive.h /$HOME/include
#
gcc -c -g collatz_recursive.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling collatz_recursive.c"
  exit
fi
rm compiler.txt
#
mv collatz_recursive.o ~/libc/$ARCH/collatz_recursive.o
#
echo "Library installed as ~/libc/$ARCH/collatz_recursive.o"
