#!/bin/bash
#
gcc -c -g -I/$HOME/include index_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling index_prb.c"
  exit
fi
rm compiler.txt
#
gcc index_prb.o /$HOME/libc/$ARCH/index.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading index_prb.o."
  exit
fi
#
rm index_prb.o
#
mv a.out index_prb
./index_prb > index_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running index_prb."
  exit
fi
rm index_prb
#
echo "Program output written to index_prb_output.txt"
