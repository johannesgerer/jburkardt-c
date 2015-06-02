#!/bin/bash
#
cp wathen.h /$HOME/include
#
gcc -c -I/$HOME/include wathen.c
if [ $? -ne 0 ]; then
  echo "Errors compiling wathen.c"
  exit
fi
#
mv wathen.o ~/libc/$ARCH/wathen.o
#
echo "Library installed as ~/libc/$ARCH/wathen.o"
