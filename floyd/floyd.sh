#!/bin/bash
#
cp floyd.h /$HOME/include
#
gcc -c -I /$HOME/include floyd.c
if [ $? -ne 0 ]; then
  echo "Errors compiling floyd.c."
  exit
fi
#
mv floyd.o ~/libc/$ARCH/floyd.o
#
echo "Library installed as ~/libc/$ARCH/floyd.o"
