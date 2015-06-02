#!/bin/bash
#
gcc -c -I/$HOME/include cube_exactness_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling cube_exactness_prb.c"
  exit
fi
#
gcc -o cube_exactness_prb cube_exactness_prb.o /$HOME/libc/$ARCH/cube_exactness.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cube_exactness_prb.o."
  exit
fi
#
rm cube_exactness_prb.o
#
./cube_exactness_prb > cube_exactness_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cube_exactness_prb."
  exit
fi
rm cube_exactness_prb
#
echo "Program output written to cube_exactness_prb_output.txt"
