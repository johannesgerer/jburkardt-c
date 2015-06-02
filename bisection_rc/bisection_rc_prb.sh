#!/bin/bash
#
gcc -c -I/$HOME/include bisection_rc_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling bisection_rc_prb.c"
  exit
fi
#
gcc bisection_rc_prb.o /$HOME/libc/$ARCH/bisection_rc.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading bisection_rc_prb.o."
  exit
fi
#
rm bisection_rc_prb.o
#
mv a.out bisection_rc_prb
./bisection_rc_prb > bisection_rc_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running bisection_rc_prb."
  exit
fi
rm bisection_rc_prb
#
echo "Program output written to bisection_rc_prb_output.txt"
