#!/bin/bash
#
cp subpak.h /$HOME/include
#
gcc -c -g subpak.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling subpak.c."
  exit
fi
rm compiler.txt
#
mv subpak.o ~/libc/$ARCH/subpak.o
#
echo "Library installed as ~/libc/$ARCH/subpak.o"
