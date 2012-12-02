#!/bin/bash
#
gcc -c -g sandia_rules_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sandia_rules_prb.c."
  exit
fi
rm compiler.txt
#
gcc sandia_rules_prb.o /$HOME/libc/$ARCH/sandia_rules.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sandia_rules_prb.o."
  exit
fi
#
rm sandia_rules_prb.o
#
mv a.out sandia_rules_prb
./sandia_rules_prb > sandia_rules_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sandia_rules_prb."
  exit
fi
rm sandia_rules_prb
#
echo "Program output written to sandia_rules_prb_output.txt"
