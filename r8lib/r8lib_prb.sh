#!/bin/bash
#
gcc -c -I/$HOME/include r8lib_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling r8lib_prb.c"
  exit
fi
#
gcc -o r8lib_prb r8lib_prb.o /$HOME/libc/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading r8lib_prb.o."
  exit
fi
#
rm r8lib_prb.o
#
./r8lib_prb > r8lib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running r8lib_prb."
  exit
fi
rm r8lib_prb
#
echo "Program output written to r8lib_prb_output.txt"
