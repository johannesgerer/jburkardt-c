#!/bin/bash
#
gcc -c -g i4lib_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling i4lib_prb.c."
  exit
fi
rm compiler.txt
#
gcc i4lib_prb.o /$HOME/libc/$ARCH/i4lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading i4lib_prb.o."
  exit
fi
#
rm i4lib_prb.o
#
mv a.out i4lib_prb
./i4lib_prb > i4lib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running i4lib_prb."
  exit
fi
rm i4lib_prb
#
echo "Program output written to i4lib_prb_output.txt"
