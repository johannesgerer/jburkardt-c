#!/bin/bash
#
cp beta_nc.h /$HOME/include
#
gcc -c -g beta_nc.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling beta_nc.c."
  exit
fi
rm compiler.txt
#
mv beta_nc.o ~/libc/$ARCH/beta_nc.o
#
echo "Library installed as ~/libc/$ARCH/beta_nc.o"
