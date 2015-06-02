#!/bin/bash
#
cp zero_rc.h /$HOME/include
#
gcc -c -I/$HOME/include zero_rc.c
if [ $? -ne 0 ]; then
  echo "Errors compiling zero_rc.c"
  exit
fi
#
mv zero_rc.o ~/libc/$ARCH/zero_rc.o
#
echo "Library installed as ~/libc/$ARCH/zero_rc.o"
