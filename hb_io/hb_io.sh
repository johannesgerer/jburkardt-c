#!/bin/bash
#
cp hb_io.h /$HOME/include
#
gcc -c -g -I /$HOME/include hb_io.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hb_io.c"
  exit
fi
rm compiler.txt
#
mv hb_io.o ~/libc/$ARCH/hb_io.o
#
echo "Library installed as ~/libc/$ARCH/hb_io.o"
