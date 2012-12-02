#!/bin/bash
#
gcc -c -I/usr/local/dislin dislin_ex11.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dislin_ex11.c"
  exit
fi
rm compiler.txt
#
gcc dislin_ex11.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dislin_ex11.o."
  exit
fi
#
rm dislin_ex11.o
#
mv a.out dislin_ex11
./dislin_ex11 > dislin_ex11_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dislin_ex11."
  exit
fi
rm dislin_ex11
#
echo "Program output written to dislin_ex11_output.txt"
