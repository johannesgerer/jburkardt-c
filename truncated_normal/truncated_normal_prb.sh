#!/bin/bash
#
gcc -c -I/$HOME/include truncated_normal_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling truncated_normal_prb.c"
  exit
fi
#
gcc truncated_normal_prb.o /$HOME/libc/$ARCH/truncated_normal.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading truncated_normal_prb.o."
  exit
fi
#
rm truncated_normal_prb.o
#
mv a.out truncated_normal_prb
./truncated_normal_prb > truncated_normal_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running truncated_normal_prb."
  exit
fi
rm truncated_normal_prb
#
echo "Program output written to truncated_normal_prb_output.txt"
