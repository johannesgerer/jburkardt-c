#!/bin/bash
#
gcc -c -I/$HOME/include walsh_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling walsh_prb.c."
  exit
fi
#
gcc walsh_prb.o /$HOME/libc/$ARCH/walsh.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading walsh_prb.o."
  exit
fi
#
rm walsh_prb.o
#
mv a.out walsh_prb
./walsh_prb > walsh_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running walsh_prb."
  exit
fi
rm walsh_prb
#
echo "Program output written to walsh_prb_output.txt"
