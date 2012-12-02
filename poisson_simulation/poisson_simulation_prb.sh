#!/bin/bash
#
gcc -c -g -I/$HOME/include poisson_simulation_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling poisson_simulation_prb.c"
  exit
fi
rm compiler.txt
#
gcc poisson_simulation_prb.o /$HOME/libc/$ARCH/poisson_simulation.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading poisson_simulation_prb.o"
  exit
fi
#
rm poisson_simulation_prb.o
#
mv a.out poisson_simulation_prb
./poisson_simulation_prb > poisson_simulation_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running poisson_simulation_prb."
  exit
fi
rm poisson_simulation_prb
#
echo "Program output written to poisson_simulation_prb_output.txt"
