#!/bin/bash
#
gcc -c -I/$HOME/include line_fekete_rule_prb.c
if [ $? -ne 0 ]; then
  echo "Errors compiling line_fekete_rule_prb.c"
  exit
fi
#
gcc line_fekete_rule_prb.o /$HOME/libc/$ARCH/line_fekete_rule.o /$HOME/libc/$ARCH/qr_solve.o /$HOME/libc/$ARCH/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading line_fekete_rule_prb.o."
  exit
fi
#
rm line_fekete_rule_prb.o
#
mv a.out line_fekete_rule_prb
./line_fekete_rule_prb > line_fekete_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running line_fekete_rule_prb."
  exit
fi
rm line_fekete_rule_prb
#
echo "Program output written to line_fekete_rule_prb_output.txt"
