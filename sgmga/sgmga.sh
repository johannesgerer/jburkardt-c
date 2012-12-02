#!/bin/bash
#
cp sgmga.h /$HOME/include
#
gcc -c -g -I /$HOME/include sgmga.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sgmga.c."
  exit
fi
rm compiler.txt
#
mv sgmga.o ~/libc/$ARCH/sgmga.o
#
echo "Library installed as ~/libc/$ARCH/sgmga.o"
