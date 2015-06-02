#!/bin/bash
#
gcc -c -I/$HOME/include square_arbq_rule_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling square_arbq_rule_prb.c"
  exit
fi
#
gcc -o square_arbq_rule_prb square_arbq_rule_prb.o /$HOME/libc/$ARCH/square_arbq_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading square_arbq_rule_prb.o."
  exit
fi
#
rm square_arbq_rule_prb.o
#
./square_arbq_rule_prb > square_arbq_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running square_arbq_rule_prb."
  exit
fi
rm square_arbq_rule_prb
#
echo "Program output written to square_arbq_rule_prb_output.txt"
