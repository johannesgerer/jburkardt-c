#!/bin/bash
#
gcc -c pointers.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pointers.c."
  exit
fi
rm compiler.txt
#
gcc pointers.o
if [ $? -ne 0 ]; then
  echo "Errors linking pointers.o."
  exit
fi
#
rm pointers.o
#
mv a.out pointers
./pointers > pointers_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pointers."
  exit
fi
rm pointers
#
echo "Program output written to pointers_output.txt"
