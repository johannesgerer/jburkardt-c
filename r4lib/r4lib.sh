#!/bin/bash
#
cp r4lib.h /$HOME/include
#
gcc -c -I /$HOME/include r4lib.c
if [ $? -ne 0 ]; then
  echo "Errors compiling r4lib.c"
  exit
fi
#
mv r4lib.o ~/libc/$ARCH/r4lib.o
#
echo "Library installed as ~/libc/$ARCH/r4lib.o"
