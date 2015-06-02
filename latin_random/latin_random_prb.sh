#!/bin/bash
#
gcc -c -I/$HOME/include latin_random_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling latin_random_prb.c"
  exit
fi
#
gcc latin_random_prb.o /$HOME/libc/$ARCH/latin_random.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading latin_random_prb.o."
  exit
fi
#
rm latin_random_prb.o
#
mv a.out latin_random_prb
./latin_random_prb > latin_random_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running latin_random_prb."
  exit
fi
rm latin_random_prb
#
echo "Program output written to latin_random_prb_output.txt"
