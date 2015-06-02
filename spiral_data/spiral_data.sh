#!/bin/bash
#
cp spiral_data.h /$HOME/include
#
gcc -c -I/$HOME/include spiral_data.c
if [ $? -ne 0 ]; then
  echo "Errors compiling spiral_data.c"
  exit
fi
#
mv spiral_data.o ~/libc/$ARCH/spiral_data.o
#
echo "Library installed as ~/libc/$ARCH/spiral_data.o"
