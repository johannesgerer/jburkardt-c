#!/bin/bash
#
gcc -c -g quadrule_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quadrule_prb.c."
  exit
fi
rm compiler.txt
#
gcc quadrule_prb.o /$HOME/libc/$ARCH/quadrule.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quadrule_prb.o."
  exit
fi
#
rm quadrule_prb.o
#
mv a.out quadrule_prb
./quadrule_prb > quadrule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running quadrule_prb."
  exit
fi
rm quadrule_prb
#
echo "Program output written to quadrule_prb_output.txt"
