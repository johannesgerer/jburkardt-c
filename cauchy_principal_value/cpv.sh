#!/bin/bash
#
cp cpv.h /$HOME/include
#
gcc -c -I/$HOME/include cpv.c
if [ $? -ne 0 ]; then
  echo "Errors compiling cpv.c"
  exit
fi
#
mv cpv.o ~/libc/$ARCH/cpv.o
#
echo "Library installed as ~/libc/$ARCH/cpv.o"
