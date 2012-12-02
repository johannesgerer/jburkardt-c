#!/bin/bash
#
cp task_division.h /$HOME/include
#
gcc -c -g -I /$HOME/include task_division.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling task_division.c."
  exit
fi
rm compiler.txt
#
mv task_division.o ~/libc/$ARCH/task_division.o
#
echo "Library installed as ~/libc/$ARCH/task_division.o"
