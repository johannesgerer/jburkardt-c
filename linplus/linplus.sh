#!/bin/bash
#
cp linplus.h /$HOME/include
#
gcc -c linplus.c
if [ $? -ne 0 ]; then
  echo "Errors compiling linplus.c."
  exit
fi
#
mv linplus.o ~/libc/$ARCH/linplus.o
#
echo "Library installed as ~/libc/$ARCH/linplus.o"
