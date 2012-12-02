#!/bin/bash
#
cp compass_search.h /$HOME/include
#
gcc -c -g -I /$HOME/include compass_search.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling compass_search.c"
  exit
fi
rm compiler.txt
#
mv compass_search.o ~/libc/$ARCH/compass_search.o
#
echo "Library installed as ~/libc/$ARCH/compass_search.o"
