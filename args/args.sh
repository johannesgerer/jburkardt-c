#!/bin/bash
#
gcc -c -g args.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling args.c"
  exit
fi
rm compiler.txt
#
gcc args.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading args.o"
  exit
fi
#
rm args.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/args
#
echo "Executable installed as ~/binc/$ARCH/args"
