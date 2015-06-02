#!/bin/bash
#
cp line_integrals.h /$HOME/include
#
gcc -c -g -I /$HOME/include line_integrals.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling line_integrals.c"
  exit
fi
rm compiler.txt
#
mv line_integrals.o ~/libc/$ARCH/line_integrals.o
#
echo "Library installed as ~/libc/$ARCH/line_integrals.o"
