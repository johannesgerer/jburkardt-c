#!/bin/bash
#
cp rnglib.h /$HOME/include
#
gcc -c -I/$HOME/include rnglib.c
if [ $? -ne 0 ]; then
  echo "Errors compiling rnglib.c"
  exit
fi
#
mv rnglib.o ~/libc/$ARCH/rnglib.o
#
echo "Library installed as ~/libc/$ARCH/rnglib.o"
