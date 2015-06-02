#!/bin/bash
#
gcc -c -I/$HOME/include cc_io_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling cc_io_prb.c"
  exit
fi
#
gcc -o cc_io_prb cc_io_prb.o /$HOME/libc/$ARCH/cc_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cc_io_prb.o."
  exit
fi
#
rm cc_io_prb.o
#
./cc_io_prb > cc_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cc_io_prb."
  exit
fi
rm cc_io_prb
#
echo "Program output written to cc_io_prb_output.txt"
