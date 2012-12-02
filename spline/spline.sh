#!/bin/bash
#
cp spline.h /$HOME/include
#
gcc -c -g -I /$HOME/include spline.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling spline.c"
  exit
fi
rm compiler.txt
#
mv spline.o ~/libc/$ARCH/spline.o
#
echo "Library installed as ~/libc/$ARCH/spline.o"
