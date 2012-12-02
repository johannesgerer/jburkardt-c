#!/bin/bash
#
cp weekday.h /$HOME/include
#
gcc -c -g -I /$HOME/include weekday.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling weekday.c"
  exit
fi
rm compiler.txt
#
mv weekday.o ~/libc/$ARCH/weekday.o
#
echo "Library installed as ~/libc/$ARCH/weekday.o"
