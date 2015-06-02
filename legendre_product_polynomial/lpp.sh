#!/bin/bash
#
cp lpp.h /$HOME/include
#
gcc -c -I/$HOME/include lpp.c
if [ $? -ne 0 ]; then
  echo "Errors compiling lpp.c"
  exit
fi
#
mv lpp.o ~/libc/$ARCH/lpp.o
#
echo "Library installed as ~/libc/$ARCH/lpp.o"
