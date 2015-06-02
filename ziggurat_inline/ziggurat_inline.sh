#!/bin/bash
#
cp ziggurat_inline.h /$HOME/include
#
gcc -c ziggurat_inline.c
if [ $? -ne 0 ]; then
  echo "Errors compiling ziggurat_inline.c."
  exit
fi
#
mv ziggurat_inline.o ~/libc/$ARCH/ziggurat_inline.o
#
echo "Library installed as ~/libc/$ARCH/ziggurat_inline.o"
