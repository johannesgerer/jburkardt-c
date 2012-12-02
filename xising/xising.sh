#!/bin/bash
#
gcc -c xising.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling xising.c"
  exit
fi
rm compiler.txt
#
gcc -O -L/usr/X11R6/lib xising.c -lX11 
if [ $? -ne 0 ]; then
  echo "Errors linking and loading xising.o"
  exit
fi
#
rm xising.o
mv a.out ~/binc/$ARCH/xising
#
echo "Program installed as ~/binc/$ARCH/xising"
