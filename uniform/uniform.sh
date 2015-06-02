#!/bin/bash
#
cp uniform.h /$HOME/include
#
gcc -c uniform.c
if [ $? -ne 0 ]; then
  echo "Errors compiling uniform.c."
  exit
fi
#
mv uniform.o ~/libc/$ARCH/uniform.o
#
echo "Library installed as ~/libc/$ARCH/uniform.o"
