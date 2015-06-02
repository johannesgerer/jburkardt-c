#!/bin/bash
#
gcc -c -I/$HOME/include lpp_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling lpp_prb.c"
  exit
fi
#
gcc -o lpp_prb lpp_prb.o /$HOME/libc/$ARCH/lpp.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lpp_prb.o."
  exit
fi
#
rm lpp_prb.o
#
./lpp_prb > lpp_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lpp_prb."
  exit
fi
rm lpp_prb
#
echo "Program output written to lpp_prb_output.txt"
