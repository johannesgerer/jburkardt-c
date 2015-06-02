#!/bin/bash
#
cp cc_io.h /$HOME/include
#
gcc -c -I/$HOME/include cc_io.c
if [ $? -ne 0 ]; then
  echo "Errors compiling cc_io.c"
  exit
fi
#
mv cc_io.o ~/libc/$ARCH/cc_io.o
#
echo "Library installed as ~/libc/$ARCH/cc_io.o"
