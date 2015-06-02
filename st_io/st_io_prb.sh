#!/bin/bash
#
gcc -c -I/$HOME/include st_io_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling st_io_prb.c"
  exit
fi
#
gcc st_io_prb.o /$HOME/libc/$ARCH/st_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading st_io_prb.o."
  exit
fi
#
rm st_io_prb.o
#
mv a.out st_io_prb
./st_io_prb > st_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running st_io_prb."
  exit
fi
rm st_io_prb
#
echo "Program output written to st_io_prb_output.txt"
