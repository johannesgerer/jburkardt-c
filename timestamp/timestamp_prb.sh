#!/bin/bash
#
gcc -c -g timestamp_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling timestamp_prb.c."
  exit
fi
rm compiler.txt
#
gcc timestamp_prb.o /$HOME/libc/$ARCH/timestamp.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading timestamp_prb.o."
  exit
fi
#
rm timestamp_prb.o
#
mv a.out timestamp_prb
./timestamp_prb > timestamp_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running timestamp_prb."
  exit
fi
rm timestamp_prb
#
echo "Program output written to timestamp_prb_output.txt"
