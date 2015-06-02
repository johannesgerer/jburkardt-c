#!/bin/bash
#
cp sort_rc.h /$HOME/include
#
gcc -c -I/$HOME/include sort_rc.c
if [ $? -ne 0 ]; then
  echo "Errors compiling sort_rc.c"
  exit
fi
#
mv sort_rc.o ~/libc/$ARCH/sort_rc.o
#
echo "Library installed as ~/libc/$ARCH/sort_rc.o"
