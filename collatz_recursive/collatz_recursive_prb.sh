#!/bin/bash
#
gcc -c -g collatz_recursive_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling collatz_recursive_prb.c"
  exit
fi
rm compiler.txt
#
gcc collatz_recursive_prb.o /$HOME/libc/$ARCH/collatz_recursive.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading collatz_recursive_prb.o."
  exit
fi
#
rm collatz_recursive_prb.o
#
mv a.out collatz_recursive_prb
./collatz_recursive_prb > collatz_recursive_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running collatz_recursive_prb."
  exit
fi
rm collatz_recursive_prb
#
echo "Program output written to collatz_recursive_prb_output.txt"
