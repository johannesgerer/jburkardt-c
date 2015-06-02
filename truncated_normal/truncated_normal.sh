#!/bin/bash
#
cp truncated_normal.h /$HOME/include
#
gcc -c -I/$HOME/include truncated_normal.c
if [ $? -ne 0 ]; then
  echo "Errors compiling truncated_normal.c"
  exit
fi
#
mv truncated_normal.o ~/libc/$ARCH/truncated_normal.o
#
echo "Library installed as ~/libc/$ARCH/truncated_normal.o"
