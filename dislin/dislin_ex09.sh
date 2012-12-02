#!/bin/bash
#
gcc -c -I/usr/local/dislin dislin_ex09.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dislin_ex09.c"
  exit
fi
rm compiler.txt
#
gcc dislin_ex09.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dislin_ex09.o."
  exit
fi
#
rm dislin_ex09.o
#
mv a.out dislin_ex09
./dislin_ex09 > dislin_ex09_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dislin_ex09."
  exit
fi
rm dislin_ex09
#
echo "Program output written to dislin_ex09_output.txt"
