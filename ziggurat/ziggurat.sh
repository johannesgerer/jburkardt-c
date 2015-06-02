#!/bin/bash
#
cp ziggurat.h /$HOME/include
#
gcc -c ziggurat.c
if [ $? -ne 0 ]; then
  echo "Errors compiling ziggurat.c."
  exit
fi
#
mv ziggurat.o ~/libc/$ARCH/ziggurat.o
#
echo "Library installed as ~/libc/$ARCH/ziggurat.o"
