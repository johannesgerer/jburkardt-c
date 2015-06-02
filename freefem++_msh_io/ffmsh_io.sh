#!/bin/bash
#
cp ffmsh_io.h /$HOME/include
#
gcc -c -I/$HOME/include ffmsh_io.c
if [ $? -ne 0 ]; then
  echo "Errors compiling ffmsh_io.c"
  exit
fi
#
mv ffmsh_io.o ~/libc/$ARCH/ffmsh_io.o
#
echo "Library installed as ~/libc/$ARCH/ffmsh_io.o"
