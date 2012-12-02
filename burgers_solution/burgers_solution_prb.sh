#!/bin/bash
#
gcc -c -g burgers_solution_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling burgers_solution_prb.c."
  exit
fi
rm compiler.txt
#
gcc burgers_solution_prb.o /$HOME/libc/$ARCH/burgers_solution.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading burgers_solution_prb.o."
  exit
fi
#
rm burgers_solution_prb.o
#
mv a.out burgers_solution_prb
./burgers_solution_prb > burgers_solution_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running burgers_solution_prb."
  exit
fi
rm burgers_solution_prb
#
echo "Program output written to burgers_solution_prb_output.txt"
