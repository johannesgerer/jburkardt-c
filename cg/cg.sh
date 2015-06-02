#!/bin/bash
#
cp cg.h /$HOME/include
#
gcc -c -I/$HOME/include cg.c
if [ $? -ne 0 ]; then
  echo "Errors compiling cg.c"
  exit
fi
#
mv cg.o ~/libc/$ARCH/cg.o
#
echo "Library installed as ~/libc/$ARCH/cg.o"
