#!/bin/bash
#
gcc -c -I/$HOME/include fem_io_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling fem_io_prb.c"
  exit
fi
#
gcc fem_io_prb.o /$HOME/libc/$ARCH/fem_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem_io_prb.o."
  exit
fi
#
rm fem_io_prb.o
#
mv a.out fem_io_prb
./fem_io_prb > fem_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fem_io_prb."
  exit
fi
rm fem_io_prb
#
echo "Program output written to fem_io_prb_output.txt"
