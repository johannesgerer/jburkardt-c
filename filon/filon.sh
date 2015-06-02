#!/bin/bash
#
cp filon.h /$HOME/include
#
gcc -c -I/$HOME/include filon.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling filon.c"
  exit
fi
rm compiler.txt
#
mv filon.o ~/libc/$ARCH/filon.o
#
echo "Library installed as ~/libc/$ARCH/filon.o"
