#!/bin/bash
#
cp c8lib.h /$HOME/include
#
gcc -c -I/$HOME/include c8lib.c
if [ $? -ne 0 ]; then
  echo "Errors compiling c8lib.c."
  exit
fi
#
mv c8lib.o ~/libc/$ARCH/c8lib.o
#
echo "Library installed as ~/libc/$ARCH/c8lib.o"
