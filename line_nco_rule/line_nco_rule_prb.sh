#!/bin/bash
#
gcc -c -I/$HOME/include line_nco_rule_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling line_nco_rule_prb.c"
  exit
fi
#
gcc line_nco_rule_prb.o /$HOME/libc/$ARCH/line_nco_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading line_nco_rule_prb.o."
  exit
fi
#
rm line_nco_rule_prb.o
#
mv a.out line_nco_rule_prb
./line_nco_rule_prb > line_nco_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running line_nco_rule_prb."
  exit
fi
rm line_nco_rule_prb
#
echo "Program output written to line_nco_rule_prb_output.txt"
