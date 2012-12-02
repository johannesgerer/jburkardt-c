#!/bin/bash
#
gcc -c -g -I/$HOME/include sine_transform_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sine_transform_prb.c"
  exit
fi
rm compiler.txt
#
gcc sine_transform_prb.o /$HOME/libc/$ARCH/sine_transform.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sine_transform_prb.o."
  exit
fi
#
rm sine_transform_prb.o
#
mv a.out sine_transform_prb
./sine_transform_prb > sine_transform_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sine_transform_prb."
  exit
fi
rm sine_transform_prb
#
echo "Program output written to sine_transform_prb_output.txt"
