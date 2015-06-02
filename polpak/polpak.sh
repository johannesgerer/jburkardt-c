#!/bin/bash
#
cp polpak.h /$HOME/include
#
gcc -c polpak.c
if [ $? -ne 0 ]; then
  echo "Errors compiling polpak.c."
  exit
fi
#
mv polpak.o ~/libc/$ARCH/polpak.o
#
echo "Library installed as ~/libc/$ARCH/polpak.o"
