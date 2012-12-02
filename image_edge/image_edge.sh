#!/bin/bash
#
cp image_edge.h /$HOME/include
#
gcc -c -g image_edge.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling image_edge.c."
  exit
fi
rm compiler.txt
#
mv image_edge.o ~/libc/$ARCH/image_edge.o
#
echo "Library installed as ~/libc/$ARCH/image_edge.o"
