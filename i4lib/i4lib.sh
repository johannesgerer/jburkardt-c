#!/bin/bash
#
cp i4lib.h /$HOME/include
#
gcc -c i4lib.c
if [ $? -ne 0 ]; then
  echo "Errors compiling i4lib.c."
  exit
fi
#
mv i4lib.o ~/libc/$ARCH/i4lib.o
#
echo "Library installed as ~/libc/$ARCH/i4lib.o"
