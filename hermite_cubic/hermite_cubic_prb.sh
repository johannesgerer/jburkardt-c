#!/bin/bash
#
gcc -c -g -I/$HOME/include hermite_cubic_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hermite_cubic_prb.c"
  exit
fi
rm compiler.txt
#
gcc hermite_cubic_prb.o /$HOME/libc/$ARCH/hermite_cubic.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hermite_cubic_prb.o."
  exit
fi
#
rm hermite_cubic_prb.o
#
mv a.out hermite_cubic_prb
./hermite_cubic_prb > hermite_cubic_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running hermite_cubic_prb."
  exit
fi
rm hermite_cubic_prb
#
echo "Program output written to hermite_cubic_prb_output.txt"
