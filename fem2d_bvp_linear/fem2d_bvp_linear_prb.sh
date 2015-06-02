#!/bin/bash
#
gcc -c fem2d_bvp_linear_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling fem2d_bvp_linear_prb.c"
  exit
fi
#
gcc fem2d_bvp_linear_prb.o ~/libc/$ARCH/fem2d_bvp_linear.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem2d_bvp_linear_prb.o"
  exit
fi
rm fem2d_bvp_linear_prb.o
#
mv a.out fem2d_bvp_linear_prb
./fem2d_bvp_linear_prb > fem2d_bvp_linear_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fem2d_bvp_linear_prb"
  exit
fi
rm fem2d_bvp_linear_prb
#
echo "Program output written to fem2d_bvp_linear_prb_output.txt"
