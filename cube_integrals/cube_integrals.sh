#!/bin/bash
#
cp cube_integrals.h /$HOME/include
#
gcc -c -g -I /$HOME/include cube_integrals.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cube_integrals.c"
  exit
fi
rm compiler.txt
#
mv cube_integrals.o ~/libc/$ARCH/cube_integrals.o
#
echo "Library installed as ~/libc/$ARCH/cube_integrals.o"
