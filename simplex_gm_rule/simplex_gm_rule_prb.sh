#!/bin/bash
#
gcc -c -I/$HOME/include simplex_gm_rule_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling simplex_gm_rule_prb.c"
  exit
fi
#
gcc simplex_gm_rule_prb.o /$HOME/libc/$ARCH/simplex_gm_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading simplex_gm_rule_prb.o."
  exit
fi
#
rm simplex_gm_rule_prb.o
#
mv a.out simplex_gm_rule_prb
./simplex_gm_rule_prb > simplex_gm_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running simplex_gm_rule_prb."
  exit
fi
rm simplex_gm_rule_prb
#
echo "Program output written to simplex_gm_rule_prb_output.txt"
