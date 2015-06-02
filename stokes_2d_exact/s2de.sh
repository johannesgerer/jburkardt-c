#!/bin/bash
#
cp s2de.h /$HOME/include
#
gcc -c -I/$HOME/include s2de.c
if [ $? -ne 0 ]; then
  echo "Errors compiling s2de.c"
  exit
fi
#
mv s2de.o ~/libc/$ARCH/s2de.o
#
echo "Library installed as ~/libc/$ARCH/s2de.o"
