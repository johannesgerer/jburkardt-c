#!/bin/bash
#
gcc -c -g problem3.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling problem3.c"
  exit
fi
rm compiler.txt
#
gcc problem3.o $HOME/libc/$ARCH/fd1d_heat_steady.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading problem3.o"
  exit
fi
rm problem3.o
#
mv a.out problem3
./problem3 > problem3_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running problem3"
  exit
fi
rm problem3
#
echo "Test program output written to problem3_output.txt."
