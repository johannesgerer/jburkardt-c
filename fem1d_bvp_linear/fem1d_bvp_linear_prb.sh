#!/bin/bash
#
gcc -c fem1d_bvp_linear_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_bvp_linear_prb.c"
  exit
fi
#
gcc fem1d_bvp_linear_prb.o ~/libc/$ARCH/fem1d_bvp_linear.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem1d_bvp_linear_prb.o"
  exit
fi
rm fem1d_bvp_linear_prb.o
#
mv a.out fem1d_bvp_linear_prb
./fem1d_bvp_linear_prb > fem1d_bvp_linear_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fem1d_bvp_linear_prb"
  exit
fi
rm fem1d_bvp_linear_prb
#
echo "Program output written to fem1d_bvp_linear_prb_output.txt"
