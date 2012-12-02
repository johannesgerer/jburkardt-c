#!/bin/bash
#
gcc -c -g pgma_io_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pgma_io_prb.c."
  exit
fi
rm compiler.txt
#
gcc pgma_io_prb.o /$HOME/libc/$ARCH/pgma_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pgma_io_prb.o."
  exit
fi
#
rm pgma_io_prb.o
#
mv a.out pgma_io_prb
./pgma_io_prb > pgma_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pgma_io_prb."
  exit
fi
rm pgma_io_prb
#
echo "Program output written to pgma_io_prb_output.txt"
