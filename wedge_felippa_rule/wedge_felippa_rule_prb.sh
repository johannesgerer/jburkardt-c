#!/bin/bash
#
gcc -c -I/$HOME/include wedge_felippa_rule_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling wedge_felippa_rule_prb.c"
  exit
fi
#
gcc -o wedge_felippa_rule_prb wedge_felippa_rule_prb.o /$HOME/libc/$ARCH/wedge_felippa_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading wedge_felippa_rule_prb.o."
  exit
fi
#
rm wedge_felippa_rule_prb.o
#
./wedge_felippa_rule_prb > wedge_felippa_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running wedge_felippa_rule_prb."
  exit
fi
rm wedge_felippa_rule_prb
#
echo "Program output written to wedge_felippa_rule_prb_output.txt"
