#!/bin/bash
#
cp change_making.h /$HOME/include
#
gcc -c -I/$HOME/include change_making.c
if [ $? -ne 0 ]; then
  echo "Errors compiling change_making.c"
  exit
fi
#
mv change_making.o ~/libc/$ARCH/change_making.o
#
echo "Library installed as ~/libc/$ARCH/change_making.o"
