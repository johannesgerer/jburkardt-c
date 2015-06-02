#!/bin/bash
#
gcc -c -g -I/$HOME/include random_data_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling random_data_prb.c"
  exit
fi
rm compiler.txt
#
gcc random_data_prb.o /$HOME/libc/$ARCH/random_data.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading random_data_prb.o."
  exit
fi
#
rm random_data_prb.o
#
mv a.out random_data_prb
./random_data_prb > random_data_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running random_data_prb."
  exit
fi
rm random_data_prb
#
echo "Program output written to random_data_prb_output.txt"
