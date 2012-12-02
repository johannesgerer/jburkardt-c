#!/bin/bash
#
gcc -c -g correlation_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling correlation_prb.c"
  exit
fi
rm compiler.txt
#
gcc correlation_prb.o ~/libc/$ARCH/correlation.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading correlation_prb.o"
  exit
fi
rm correlation_prb.o
#
mv a.out correlation_prb
./correlation_prb > correlation_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running correlation_prb"
  exit
fi
rm correlation_prb
#
echo "Test program output written to correlation_prb_output.txt."
