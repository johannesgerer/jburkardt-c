#!/bin/bash
#
gcc -c rbox.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rbox.c."
  exit
fi
rm compiler.txt
#
gcc rbox.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading rbox.o."
  exit
fi
#
rm rbox.o
#
chmod u+x a.out
mv a.out ~/binc/$ARCH/rbox
#
echo "Executable installed as ~/binc/$ARCH/rbox"
