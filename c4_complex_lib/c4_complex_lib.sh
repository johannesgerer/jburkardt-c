#!/bin/bash
#
cp c4_complex_lib.h /$HOME/include
#
gcc -c -g c4_complex_lib.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling c4_complex_lib.c."
  exit
fi
rm compiler.txt
#
mv c4_complex_lib.o ~/libc/$ARCH/c4_complex_lib.o
#
echo "Library installed as ~/libc/$ARCH/c4_complex_lib.o"
