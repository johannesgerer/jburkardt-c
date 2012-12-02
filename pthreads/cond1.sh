#!/bin/bash
#
gcc -c -g cond1.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cond1.c"
  exit
fi
rm compiler.txt
#
gcc cond1.o
if [ $? -ne 0 ]; then
  echo "Errors loading cond1.c"
  exit
fi
#
rm cond1.o
mv a.out cond1
cond1 > cond1_output.txt
rm cond1
#
echo "Program output written to cond1_output.txt"
