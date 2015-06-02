#!/bin/bash
#
cp hpp.h /$HOME/include
#
gcc -c -I/$HOME/include hpp.c
if [ $? -ne 0 ]; then
  echo "Errors compiling hpp.c"
  exit
fi
#
mv hpp.o ~/libc/$ARCH/hpp.o
#
echo "Library installed as ~/libc/$ARCH/hpp.o"
