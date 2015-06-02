#!/bin/bash
#
cp spline.h /$HOME/include
#
gcc -c -I /$HOME/include spline.c
if [ $? -ne 0 ]; then
  echo "Errors compiling spline.c"
  exit
fi
#
mv spline.o ~/libc/$ARCH/spline.o
#
echo "Library installed as ~/libc/$ARCH/spline.o"
