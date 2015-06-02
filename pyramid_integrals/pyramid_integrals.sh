#!/bin/bash
#
cp pyramid_integrals.h /$HOME/include
#
gcc -c -I/$HOME/include pyramid_integrals.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pyramid_integrals.c"
  exit
fi
rm compiler.txt
#
mv pyramid_integrals.o ~/libc/$ARCH/pyramid_integrals.o
#
echo "Library installed as ~/libc/$ARCH/pyramid_integrals.o"
