#!/bin/bash
#
cp ranlib.h /$HOME/include
#
gcc -c -I/$HOME/include ranlib.c
if [ $? -ne 0 ]; then
  echo "Errors compiling ranlib.c"
  exit
fi
#
mv ranlib.o ~/libc/$ARCH/ranlib.o
#
echo "Library installed as ~/libc/$ARCH/ranlib.o"
