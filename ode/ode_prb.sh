#!/bin/bash
#
gcc -c -g -I/$HOME/include ode_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ode_prb.c"
  exit
fi
rm compiler.txt
#
gcc ode_prb.o /$HOME/libc/$ARCH/ode.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ode_prb.o."
  exit
fi
#
rm ode_prb.o
#
mv a.out ode_prb
./ode_prb > ode_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ode_prb."
  exit
fi
rm ode_prb
#
echo "Program output written to ode_prb_output.txt"
