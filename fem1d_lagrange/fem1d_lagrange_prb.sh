#!/bin/bash
#
gcc -c -I/$HOME/include fem1d_lagrange_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_lagrange_prb.c"
  exit
fi
#
gcc -o fem1d_lagrange_prb fem1d_lagrange_prb.o /$HOME/libc/$ARCH/fem1d_lagrange.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem1d_lagrange_prb.o."
  exit
fi
#
rm fem1d_lagrange_prb.o
#
./fem1d_lagrange_prb > fem1d_lagrange_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fem1d_lagrange_prb."
  exit
fi
rm fem1d_lagrange_prb
#
echo "Program output written to fem1d_lagrange_prb_output.txt"
