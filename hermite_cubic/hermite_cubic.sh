#!/bin/bash
#
cp hermite_cubic.h /$HOME/include
#
gcc -c -g -I /$HOME/include hermite_cubic.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hermite_cubic.c"
  exit
fi
rm compiler.txt
#
mv hermite_cubic.o ~/libc/$ARCH/hermite_cubic.o
#
echo "Library installed as ~/libc/$ARCH/hermite_cubic.o"
