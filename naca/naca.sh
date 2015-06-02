#!/bin/bash
#
cp naca.h /$HOME/include
#
gcc -c -I/$HOME/include naca.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling naca.c"
  exit
fi
rm compiler.txt
#
mv naca.o ~/libc/$ARCH/naca.o
#
echo "Library installed as ~/libc/$ARCH/naca.o"
