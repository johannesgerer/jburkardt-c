#!/bin/bash
#
cp machar.h /$HOME/include
#
gcc -c machar.c
if [ $? -ne 0 ]; then
  echo "Errors compiling machar.c."
  exit
fi
#
mv machar.o ~/libc/$ARCH/machar.o
#
echo "Library installed as ~/libc/$ARCH/machar.o"
