#!/bin/bash
#
gcc -c -g problem4.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling problem4.c"
  exit
fi
rm compiler.txt
#
gcc problem4.o $HOME/libc/$ARCH/fd1d_heat_steady.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading problem4.o"
  exit
fi
rm problem4.o
#
mv a.out problem4
./problem4 > problem4_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running problem4"
  exit
fi
rm problem4
#
echo "Test program output written to problem4_output.txt."
