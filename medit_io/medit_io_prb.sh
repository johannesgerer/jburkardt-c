#!/bin/bash
#
gcc -c -I/$HOME/include medit_io_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling medit_io_prb.c."
  exit
fi
#
gcc medit_io_prb.o /$HOME/libc/$ARCH/medit_io.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading medit_io_prb.o."
  exit
fi
#
rm medit_io_prb.o
#
mv a.out medit_io_prb
./medit_io_prb > medit_io_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running medit_io_prb."
  exit
fi
rm medit_io_prb
#
echo "Program output written to medit_io_prb_output.txt"
