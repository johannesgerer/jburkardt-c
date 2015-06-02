#!/bin/bash
#
gcc -c wtime_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling wtime_prb.c."
  exit
fi
#
gcc wtime_prb.o /$HOME/libc/$ARCH/wtime.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading wtime_prb.o."
  exit
fi
#
rm wtime_prb.o
#
mv a.out wtime_prb
./wtime_prb > wtime_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running wtime_prb."
  exit
fi
rm wtime_prb
#
echo "Program output written to wtime_prb_output.txt"
