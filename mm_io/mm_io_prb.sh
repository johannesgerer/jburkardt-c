#!/bin/bash
#
gcc -c -g mm_io_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mm_io_prb.c."
  exit
fi
rm compiler.txt
#
gcc mm_io_prb.o /$HOME/libc/$ARCH/mm_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mm_io_prb.o."
  exit
fi
#
rm mm_io_prb.o
#
mv a.out mm_io_prb
./mm_io_prb > mm_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running mm_io_prb."
  exit
fi
rm mm_io_prb
#
echo "Program output written to mm_io_prb_output.txt"
