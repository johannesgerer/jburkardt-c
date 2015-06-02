#!/bin/bash
#
gcc -c -I/$HOME/include hpp_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling hpp_prb.c"
  exit
fi
#
gcc -o hpp_prb hpp_prb.o /$HOME/libc/$ARCH/hpp.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hpp_prb.o."
  exit
fi
#
rm hpp_prb.o
#
./hpp_prb > hpp_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running hpp_prb."
  exit
fi
rm hpp_prb
#
echo "Program output written to hpp_prb_output.txt"
