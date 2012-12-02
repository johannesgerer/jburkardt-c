#!/bin/bash
#
gcc -c -g -I/$HOME/include combination_lock_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling combination_lock_prb.c"
  exit
fi
rm compiler.txt
#
gcc combination_lock_prb.o /$HOME/libc/$ARCH/combination_lock.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading combination_lock_prb.o."
  exit
fi
#
rm combination_lock_prb.o
#
mv a.out combination_lock_prb
./combination_lock_prb > combination_lock_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running combination_lock_prb."
  exit
fi
rm combination_lock_prb
#
echo "Program output written to combination_lock_prb_output.txt"
