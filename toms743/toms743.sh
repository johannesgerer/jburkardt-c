#!/bin/bash
#
cp toms743.h /$HOME/include
#
gcc -c -I/$HOME/include toms743.c
if [ $? -ne 0 ]; then
  echo "Errors compiling toms743.c"
  exit
fi
#
mv toms743.o ~/libc/$ARCH/toms743.o
#
echo "Library installed as ~/libc/$ARCH/toms743.o"
