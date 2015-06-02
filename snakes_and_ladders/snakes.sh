#!/bin/bash
#
cp snakes.h /$HOME/include
#
gcc -c -I/$HOME/include snakes.c
if [ $? -ne 0 ]; then
  echo "Errors compiling snakes.c"
  exit
fi
#
mv snakes.o ~/libc/$ARCH/snakes.o
#
echo "Library installed as ~/libc/$ARCH/snakes.o"
