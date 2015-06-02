#!/bin/bash
#
cp st_to_cc.h /$HOME/include
#
gcc -c -I/$HOME/include st_to_cc.c
if [ $? -ne 0 ]; then
  echo "Errors compiling st_to_cc.c"
  exit
fi
#
mv st_to_cc.o ~/libc/$ARCH/st_to_cc.o
#
echo "Library installed as ~/libc/$ARCH/st_to_cc.o"
