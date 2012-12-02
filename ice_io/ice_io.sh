#!/bin/bash
#
cp ice_io.h /$HOME/include
#
gcc -c -g ice_io.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ice_io.c."
  exit
fi
rm compiler.txt
#
mv ice_io.o ~/libc/$ARCH/ice_io.o
#
echo "Library installed as ~/libc/$ARCH/ice_io.o"
