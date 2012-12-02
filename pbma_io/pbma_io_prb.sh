#!/bin/bash
#
gcc -c -g pbma_io_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pbma_io_prb.c."
  exit
fi
rm compiler.txt
#
gcc pbma_io_prb.o /$HOME/libc/$ARCH/pbma_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pbma_io_prb.o."
  exit
fi
#
rm pbma_io_prb.o
#
mv a.out pbma_io_prb
./pbma_io_prb > pbma_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pbma_io_prb."
  exit
fi
rm pbma_io_prb
#
echo "Program output written to pbma_io_prb_output.txt"
