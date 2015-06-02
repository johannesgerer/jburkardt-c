#!/bin/bash
#
gcc -c -O2 md.c
if [ $? -ne 0 ]; then
  echo "Errors compiling md.c"
  exit
fi
#
gcc md.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading md.o."
  exit
fi
#
rm md.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/md
#
echo "Executable installed as ~/binc/OSX/md"
