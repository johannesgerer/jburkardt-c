#!/bin/bash
#
cp c8_complex_lib.h /$HOME/include
#
gcc -c -g c8_complex_lib.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling c8_complex_lib.c."
  exit
fi
rm compiler.txt
#
mv c8_complex_lib.o ~/libc/$ARCH/c8_complex_lib.o
#
echo "Library installed as ~/libc/$ARCH/c8_complex_lib.o"
