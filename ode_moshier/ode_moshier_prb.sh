#!/bin/bash
#
gcc -c -g ode_moshier_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ode_moshier_prb.c"
  exit
fi
rm compiler.txt
#
gcc ode_moshier_prb.o -L/$HOME/libc/$ARCH -lode_moshier -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ode_moshier_prb.o."
  exit
fi
#
rm ode_moshier_prb.o
#
mv a.out ode_moshier_prb
./ode_moshier_prb < ode_moshier_prb_input.txt > ode_moshier_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running ode_moshier_prb."
  exit
fi
rm ode_moshier_prb
#
echo "Program output written to ode_moshier_prb_output.txt"
