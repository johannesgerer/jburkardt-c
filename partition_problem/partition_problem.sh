#!/bin/bash
#
cp partition_problem.h /$HOME/include
#
gcc -c -g -I /$HOME/include partition_problem.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling partition_problem.c"
  exit
fi
rm compiler.txt
#
mv partition_problem.o ~/libc/$ARCH/partition_problem.o
#
echo "Library installed as ~/libc/$ARCH/partition_problem.o"
