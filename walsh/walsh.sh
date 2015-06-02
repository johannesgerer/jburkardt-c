#!/bin/bash
#
cp walsh.h /$HOME/include
#
gcc -c -I /$HOME/include walsh.c
if [ $? -ne 0 ]; then
  echo "Errors compiling walsh.c."
  exit
fi
#
mv walsh.o ~/libc/$ARCH/walsh.o
#
echo "Library installed as ~/libc/$ARCH/walsh.o"
