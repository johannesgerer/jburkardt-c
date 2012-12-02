#!/bin/bash
#
cp geompack.h /$HOME/include
#
gcc -c -g -I /$HOME/include geompack.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling geompack.c"
  exit
fi
rm compiler.txt
#
mv geompack.o ~/libc/$ARCH/geompack.o
#
echo "Library installed as ~/libc/$ARCH/geompack.o"
