#!/bin/bash
#
gcc -c -g -I/$HOME/include brownian_motion_simulation_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling brownian_motion_simulation_prb.c"
  exit
fi
rm compiler.txt
#
gcc brownian_motion_simulation_prb.o /$HOME/libc/$ARCH/brownian_motion_simulation.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading brownian_motion_simulation_prb.o."
  exit
fi
#
rm brownian_motion_simulation_prb.o
#
mv a.out brownian_motion_simulation_prb
./brownian_motion_simulation_prb > brownian_motion_simulation_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running brownian_motion_simulation_prb."
  exit
fi
rm brownian_motion_simulation_prb
#
echo "Program output written to brownian_motion_simulation_prb_output.txt"
