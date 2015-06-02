#!/bin/bash
#
gcc -c -I/$HOME/include pyramid_integrals_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling pyramid_integrals_prb.c"
  exit
fi
rm compiler.txt
#
gcc pyramid_integrals_prb.o /$HOME/libc/$ARCH/pyramid_integrals.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pyramid_integrals_prb.o."
  exit
fi
#
rm pyramid_integrals_prb.o
#
mv a.out pyramid_integrals_prb
./pyramid_integrals_prb > pyramid_integrals_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running pyramid_integrals_prb."
  exit
fi
rm pyramid_integrals_prb
#
echo "Program output written to pyramid_integrals_prb_output.txt"
