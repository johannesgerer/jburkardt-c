#!/bin/bash
#
cp triangle_io.h /$HOME/include
#
gcc -c triangle_io.c
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_io.c."
  exit
fi
#
mv triangle_io.o ~/libc/$ARCH/triangle_io.o
#
echo "Library installed as ~/libc/$ARCH/triangle_io.o"
