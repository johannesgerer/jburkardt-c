#!/bin/bash
#
gcc -c -I/$HOME/include cg_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling cg_prb.c"
  exit
fi
#
gcc -o cg_prb cg_prb.o /$HOME/libc/$ARCH/cg.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cg_prb.o."
  exit
fi
#
rm cg_prb.o
#
./cg_prb > cg_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cg_prb."
  exit
fi
rm cg_prb
#
echo "Program output written to cg_prb_output.txt"
