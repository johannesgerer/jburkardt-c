#!/bin/bash
#
cp pink_noise.h /$HOME/include
#
gcc -c -g pink_noise.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pink_noise.c."
  exit
fi
rm compiler.txt
#
mv pink_noise.o ~/libc/$ARCH/pink_noise.o
#
echo "Library installed as ~/libc/$ARCH/pink_noise.o"
