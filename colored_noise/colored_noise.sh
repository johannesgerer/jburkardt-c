#!/bin/bash
#
cp colored_noise.h /$HOME/include
#
gcc -c -g colored_noise.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling colored_noise.c."
  exit
fi
rm compiler.txt
#
mv colored_noise.o ~/libc/$ARCH/colored_noise.o
#
echo "Library installed as ~/libc/$ARCH/colored_noise.o"
