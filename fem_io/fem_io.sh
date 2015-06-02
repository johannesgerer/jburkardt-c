#!/bin/bash
#
cp fem_io.h /$HOME/include
#
gcc -c -I /$HOME/include fem_io.c
if [ $? -ne 0 ]; then
  echo "Errors compiling fem_io.c"
  exit
fi
#
mv fem_io.o ~/libc/$ARCH/fem_io.o
#
echo "Library installed as ~/libc/$ARCH/fem_io.o"
