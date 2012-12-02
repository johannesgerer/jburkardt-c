#!/bin/bash
#
cp doomsday.h /$HOME/include
#
gcc -c -g -I /$HOME/include doomsday.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling doomsday.c"
  exit
fi
rm compiler.txt
#
mv doomsday.o ~/libc/$ARCH/doomsday.o
#
echo "Library installed as ~/libc/$ARCH/doomsday.o"
