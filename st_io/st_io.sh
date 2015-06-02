#!/bin/bash
#
cp st_io.h /$HOME/include
#
gcc -c -I/$HOME/include st_io.c
if [ $? -ne 0 ]; then
  echo "Errors compiling st_io.c"
  exit
fi
#
mv st_io.o ~/libc/$ARCH/st_io.o
#
echo "Library installed as ~/libc/$ARCH/st_io.o"
