#!/bin/bash
#
gcc -c -I/$HOME/include wathen_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling wathen_prb.c"
  exit
fi
#
gcc wathen_prb.o /$HOME/libc/$ARCH/wathen.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading wathen_prb.o."
  exit
fi
#
rm wathen_prb.o
#
mv a.out wathen_prb
./wathen_prb > wathen_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running wathen_prb."
  exit
fi
rm wathen_prb
#
echo "Program output written to wathen_prb_output.txt"
