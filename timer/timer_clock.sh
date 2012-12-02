#!/bin/bash
#
gcc -c timer_clock.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling timer_clock.c."
  exit
fi
rm compiler.txt
#
gcc timer_clock.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading timer_clock.o."
  exit
fi
#
rm timer_clock.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/timer_clock
#
echo "Executable installed as ~/binc/$ARCH/timer_clock"
