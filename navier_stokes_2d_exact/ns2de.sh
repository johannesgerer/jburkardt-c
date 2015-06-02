#!/bin/bash
#
cp ns2de.h /$HOME/include
#
gcc -c -I/$HOME/include ns2de.c
if [ $? -ne 0 ]; then
  echo "Errors compiling ns2de.c"
  exit
fi
#
mv ns2de.o ~/libc/$ARCH/ns2de.o
#
echo "Library installed as ~/libc/$ARCH/ns2de.o"
