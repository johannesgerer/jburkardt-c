#!/bin/bash
#
gcc -c -I/$HOME/include ns2de_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling ns2de_prb.c"
  exit
fi
#
gcc -o ns2de_prb ns2de_prb.o /$HOME/libc/$ARCH/ns2de.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ns2de_prb.o."
  exit
fi
#
rm ns2de_prb.o
#
./ns2de_prb > ns2de_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ns2de_prb."
  exit
fi
rm ns2de_prb
#
echo "Program output written to ns2de_prb_output.txt"
