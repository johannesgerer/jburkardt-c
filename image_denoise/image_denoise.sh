#!/bin/bash
#
cp image_denoise.h /$HOME/include
#
gcc -c -g -I /$HOME/include image_denoise.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling image_denoise.c"
  exit
fi
rm compiler.txt
#
mv image_denoise.o ~/libc/$ARCH/image_denoise.o
#
echo "Library installed as ~/libc/$ARCH/image_denoise.o"
