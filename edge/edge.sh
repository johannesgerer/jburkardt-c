#!/bin/bash
#
cp edge.h /$HOME/include
#
gcc -c -I/$HOME/include edge.c
if [ $? -ne 0 ]; then
  echo "Errors compiling edge.c"
  exit
fi
#
mv edge.o ~/libc/$ARCH/edge.o
#
echo "Library installed as ~/libc/$ARCH/edge.o"
