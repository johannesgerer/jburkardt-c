#!/bin/bash
#
cp circle_segment.h /$HOME/include
#
gcc -c -g -I/$HOME/include circle_segment.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling circle_segment.c"
  exit
fi
rm compiler.txt
#
mv circle_segment.o ~/libc/$ARCH/circle_segment.o
#
echo "Library installed as ~/libc/$ARCH/circle_segment.o"
