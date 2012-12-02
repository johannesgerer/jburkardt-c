#!/bin/bash
#
gcc -c -g  -I/$HOME/include simple_rkf45.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling simple_rkf45.c"
  exit
fi
rm compiler.txt
#
gcc simple_rkf45.o /$HOME/libc/$ARCH/rkf45.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading simple_rkf45.o"
  exit
fi
rm simple_rkf45.o
#
mv a.out simple_rkf45
./simple_rkf45 > simple_rkf45_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running simple_rkf45"
  exit
fi
rm simple_rkf45
#
echo "Test program output written to simple_rkf45_output.txt."
