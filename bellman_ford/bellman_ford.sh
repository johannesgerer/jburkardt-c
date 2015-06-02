#!/bin/bash
#
cp bellman_ford.h /$HOME/include
#
gcc -c -I /$HOME/include bellman_ford.c
if [ $? -ne 0 ]; then
  echo "Errors compiling bellman_ford.c"
  exit
fi
#
mv bellman_ford.o ~/libc/$ARCH/bellman_ford.o
#
echo "Library installed as ~/libc/$ARCH/bellman_ford.o"
