#!/bin/bash
#
cp medit_io.h /$HOME/include
#
gcc -c -I /$HOME/include medit_io.c
if [ $? -ne 0 ]; then
  echo "Errors compiling medit_io.c."
  exit
fi
#
mv medit_io.o ~/libc/$ARCH/medit_io.o
#
echo "Library installed as ~/libc/$ARCH/medit_io.o"
