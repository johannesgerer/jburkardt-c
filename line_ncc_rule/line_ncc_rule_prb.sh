#!/bin/bash
#
gcc -c -g -I/$HOME/include line_ncc_rule_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling line_ncc_rule_prb.c"
  exit
fi
rm compiler.txt
#
gcc line_ncc_rule_prb.o /$HOME/libc/$ARCH/line_ncc_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading line_ncc_rule_prb.o."
  exit
fi
#
rm line_ncc_rule_prb.o
#
mv a.out line_ncc_rule_prb
./line_ncc_rule_prb > line_ncc_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running line_ncc_rule_prb."
  exit
fi
rm line_ncc_rule_prb
#
echo "Program output written to line_ncc_rule_prb_output.txt"
