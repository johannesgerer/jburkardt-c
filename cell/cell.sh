#!/bin/bash
#
cp cell.h /$HOME/include
#
gcc -c -g -I /$HOME/include cell.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cell.c"
  exit
fi
rm compiler.txt
#
mv cell.o ~/libc/$ARCH/cell.o
#
echo "Library installed as ~/libc/$ARCH/cell.o"
