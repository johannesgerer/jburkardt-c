#!/bin/bash
#
gcc -c -g -I /$HOME/include rbf_interp_nd_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling rbf_interp_nd_prb.c"
  exit
fi
rm compiler.txt
#
gcc rbf_interp_nd_prb.o /$HOME/libc/$ARCH/rbf_interp_nd.o /$HOME/libc/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading rbf_interp_nd_prb.o."
  exit
fi
#
rm rbf_interp_nd_prb.o
#
mv a.out rbf_interp_nd_prb
./rbf_interp_nd_prb > rbf_interp_nd_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running rbf_interp_nd_prb."
  exit
fi
rm rbf_interp_nd_prb
#
echo "Program output written to rbf_interp_nd_prb_output.txt"
