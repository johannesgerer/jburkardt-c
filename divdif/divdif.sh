#!/bin/bash
#
cp divdif.h /$HOME/include
#
gcc -c -I /$HOME/include divdif.c
if [ $? -ne 0 ]; then
  echo "Errors compiling divdif.c."
  exit
fi
#
mv divdif.o ~/libc/$ARCH/divdif.o
#
echo "Library installed as ~/libc/$ARCH/divdif.o"
