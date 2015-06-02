#!/bin/bash
#
gcc -c -g -I/$HOME/include high_card_simulation_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling high_card_simulation_prb.c"
  exit
fi
rm compiler.txt
#
gcc high_card_simulation_prb.o /$HOME/libc/$ARCH/high_card_simulation.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading high_card_simulation_prb.o."
  exit
fi
#
rm high_card_simulation_prb.o
#
mv a.out high_card_simulation_prb
./high_card_simulation_prb > high_card_simulation_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running high_card_simulation_prb."
  exit
fi
rm high_card_simulation_prb
#
echo "Program output written to high_card_simulation_prb_output.txt"
