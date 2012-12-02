#!/bin/bash
#
gcc -c -g -I/$HOME/include weekday_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling weekday_prb.c"
  exit
fi
rm compiler.txt
#
gcc weekday_prb.o /$HOME/libc/$ARCH/weekday.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading weekday_prb.o."
  exit
fi
#
rm weekday_prb.o
#
mv a.out weekday_prb
./weekday_prb > weekday_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running weekday_prb."
  exit
fi
rm weekday_prb
#
echo "Program output written to weekday_prb_output.txt"
