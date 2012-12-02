#!/bin/bash
#
cp burgers_solution.h /$HOME/include
#
gcc -c -g burgers_solution.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling burgers_solution.c."
  exit
fi
rm compiler.txt
#
mv burgers_solution.o ~/libc/$ARCH/burgers_solution.o
#
echo "Library installed as ~/libc/$ARCH/burgers_solution.o"
