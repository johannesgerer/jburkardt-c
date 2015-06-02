#!/bin/bash
#
gcc -c -I/$HOME/include sort_rc_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling sort_rc_prb.c"
  exit
fi
#
gcc -o sort_rc_prb sort_rc_prb.o /$HOME/libc/$ARCH/sort_rc.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sort_rc_prb.o."
  exit
fi
#
rm sort_rc_prb.o
#
./sort_rc_prb > sort_rc_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sort_rc_prb."
  exit
fi
rm sort_rc_prb
#
echo "Program output written to sort_rc_prb_output.txt"
