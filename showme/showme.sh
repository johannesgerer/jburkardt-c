#!/bin/bash
#
gcc -c showme.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling showme.c"
  exit
fi
rm compiler.txt
#
gcc showme.o -L/usr/X11R6/lib -lX11
if [ $? -ne 0 ]; then
  echo "Errors loading showme.o"
  exit
fi
#
rm showme.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/showme
#
echo "Executable installed as ~/binc/$ARCH/showme"
