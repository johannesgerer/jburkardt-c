#!/bin/bash
#
cp padua.h /$HOME/include
#
gcc -c -I/$HOME/include padua.c
if [ $? -ne 0 ]; then
  echo "Errors compiling padua.c"
  exit
fi
#
mv padua.o ~/libc/$ARCH/padua.o
#
echo "Library installed as ~/libc/$ARCH/padua.o"
