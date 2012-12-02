#!/bin/bash
#
gcc -c -g c4_complex_lib_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling c4_complex_lib_prb.c."
  exit
fi
rm compiler.txt
#
gcc c4_complex_lib_prb.o /$HOME/libc/$ARCH/c4_complex_lib.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading c4_complex_lib_prb.o."
  exit
fi
#
rm c4_complex_lib_prb.o
#
mv a.out c4_complex_lib_prb
./c4_complex_lib_prb > c4_complex_lib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running c4_complex_lib_prb."
  exit
fi
rm c4_complex_lib_prb
#
echo "Program output written to c4_complex_lib_prb_output.txt"
