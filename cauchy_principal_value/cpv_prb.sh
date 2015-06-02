#!/bin/bash
#
gcc -c -I/$HOME/include cpv_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling cpv_prb.c"
  exit
fi
#
gcc -o cpv_prb cpv_prb.o /$HOME/libc/$ARCH/cpv.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cpv_prb.o."
  exit
fi
#
rm cpv_prb.o
#
./cpv_prb > cpv_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cpv_prb."
  exit
fi
rm cpv_prb
#
echo "Program output written to cpv_prb_output.txt"
