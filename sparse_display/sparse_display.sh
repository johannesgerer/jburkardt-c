#!/bin/bash
#
cp sparse_display.h /$HOME/include
#
gcc -c -I/$HOME/include sparse_display.c
if [ $? -ne 0 ]; then
  echo "Errors compiling sparse_display.c"
  exit
fi
#
mv sparse_display.o ~/libc/$ARCH/sparse_display.o
#
echo "Library installed as ~/libc/$ARCH/sparse_display.o"
