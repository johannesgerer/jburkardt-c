#!/bin/bash
#
gcc -c -I/$HOME/include st_to_cc_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling st_to_cc_prb.c"
  exit
fi
#
gcc -o st_to_cc_prb st_to_cc_prb.o /$HOME/libc/$ARCH/st_to_cc.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading st_to_cc_prb.o."
  exit
fi
#
rm st_to_cc_prb.o
#
./st_to_cc_prb > st_to_cc_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running st_to_cc_prb."
  exit
fi
rm st_to_cc_prb
#
echo "Program output written to st_to_cc_prb_output.txt"
