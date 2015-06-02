#!/bin/bash
#
gcc -c -I/$HOME/include snakes_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling snakes_prb.c"
  exit
fi
#
gcc -o snakes_prb snakes_prb.o /$HOME/libc/$ARCH/snakes.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading snakes_prb.o."
  exit
fi
#
rm snakes_prb.o
#
./snakes_prb > snakes_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running snakes_prb."
  exit
fi
rm snakes_prb
#
echo "Program output written to snakes_prb_output.txt"
