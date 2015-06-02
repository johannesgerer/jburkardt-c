#!/bin/bash
#
cp filum.h /$HOME/include
#
gcc -c filum.c
if [ $? -ne 0 ]; then
  echo "Errors compiling filum.c."
  exit
fi
#
mv filum.o ~/libc/$ARCH/filum.o
#
echo "Library installed as ~/libc/$ARCH/filum.o"
