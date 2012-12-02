#!/bin/bash
#
gcc -c xwaves.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling xwaves.c"
  exit
fi
rm compiler.txt
#
gcc -O -L/usr/X11R6/lib xwaves.c -lX11 
if [ $? -ne 0 ]; then
  echo "Errors linking and loading xwaves.o"
  exit
fi
#
rm xwaves.o
mv a.out ~/binc/$ARCH/xwaves
#
echo "Program installed as ~/binc/$ARCH/xwaves"
