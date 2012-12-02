#!/bin/bash
#
cp point_merge.h /$HOME/include
#
gcc -c -g point_merge.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling point_merge.c."
  exit
fi
rm compiler.txt
#
mv point_merge.o ~/libc/$ARCH/point_merge.o
#
echo "Library installed as ~/libc/$ARCH/point_merge.o"
