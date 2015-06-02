#!/bin/bash
#
gcc -c -g -I/$HOME/include c4lib_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling c4lib_prb.c."
  exit
fi
rm compiler.txt
#
gcc c4lib_prb.o /$HOME/libc/$ARCH/c4lib.o /$HOME/libc/$ARCH/r4lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading c4lib_prb.o."
  exit
fi
#
rm c4lib_prb.o
#
mv a.out c4lib_prb
./c4lib_prb > c4lib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running c4lib_prb."
  exit
fi
rm c4lib_prb
#
echo "Program output written to c4lib_prb_output.txt"
