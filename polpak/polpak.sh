#!/bin/bash
#
cp polpak.h /$HOME/include
#
gcc -c -g polpak.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling polpak.c."
  exit
fi
rm compiler.txt
#
mv polpak.o ~/libc/$ARCH/polpak.o
#
echo "Library installed as ~/libc/$ARCH/polpak.o"
