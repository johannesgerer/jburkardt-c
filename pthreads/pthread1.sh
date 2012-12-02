#!/bin/bash
#
gcc -c -g pthread1.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pthread1.c"
  exit
fi
rm compiler.txt
#
gcc pthread1.o
if [ $? -ne 0 ]; then
  echo "Errors loading pthread1.c"
  exit
fi
#
rm pthread1.o
mv a.out pthread1
pthread1 > pthread1_output.txt
rm pthread1
#
echo "Program output written to pthread1_output.txt"
