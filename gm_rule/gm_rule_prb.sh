#!/bin/bash
#
gcc -c -g -I/$HOME/include gm_rule_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling gm_rule_prb.c"
  exit
fi
rm compiler.txt
#
gcc gm_rule_prb.o /$HOME/libc/$ARCH/gm_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading gm_rule_prb.o."
  exit
fi
#
rm gm_rule_prb.o
#
mv a.out gm_rule_prb
./gm_rule_prb > gm_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running gm_rule_prb."
  exit
fi
rm gm_rule_prb
#
echo "Program output written to gm_rule_prb_output.txt"
