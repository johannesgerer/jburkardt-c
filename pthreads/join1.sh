#!/bin/bash
#
gcc -c -g join1.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling join1.c"
  exit
fi
rm compiler.txt
#
gcc join1.o
if [ $? -ne 0 ]; then
  echo "Errors loading join1.c"
  exit
fi
#
rm join1.o
mv a.out join1
join1 > join1_output.txt
rm join1
#
echo "Program output written to join1_output.txt"
