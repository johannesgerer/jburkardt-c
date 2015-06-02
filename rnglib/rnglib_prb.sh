#!/bin/bash
#
gcc -c -I/$HOME/include rnglib_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling rnglib_prb.c"
  exit
fi
#
gcc rnglib_prb.o /$HOME/libc/$ARCH/rnglib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading rnglib_prb.o."
  exit
fi
#
rm rnglib_prb.o
#
mv a.out rnglib_prb
./rnglib_prb > rnglib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running rnglib_prb."
  exit
fi
rm rnglib_prb
#
echo "Program output written to rnglib_prb_output.txt"
