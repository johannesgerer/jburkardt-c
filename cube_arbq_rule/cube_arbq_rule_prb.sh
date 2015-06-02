#!/bin/bash
#
gcc -c -I/$HOME/include cube_arbq_rule_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling cube_arbq_rule_prb.c"
  exit
fi
#
gcc -o cube_arbq_rule_prb cube_arbq_rule_prb.o /$HOME/libc/$ARCH/cube_arbq_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cube_arbq_rule_prb.o."
  exit
fi
#
rm cube_arbq_rule_prb.o
#
./cube_arbq_rule_prb > cube_arbq_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cube_arbq_rule_prb."
  exit
fi
rm cube_arbq_rule_prb
#
echo "Program output written to cube_arbq_rule_prb_output.txt"
