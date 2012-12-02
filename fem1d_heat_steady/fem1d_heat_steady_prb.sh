#!/bin/bash
#
gcc -c -g fem1d_heat_steady_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fem1d_heat_steady_prb.c"
  exit
fi
rm compiler.txt
#
gcc fem1d_heat_steady_prb.o ~/libc/$ARCH/fem1d_heat_steady.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem1d_heat_steady_prb.o"
  exit
fi
rm fem1d_heat_steady_prb.o
#
mv a.out fem1d_heat_steady_prb
./fem1d_heat_steady_prb > fem1d_heat_steady_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fem1d_heat_steady_prb"
  exit
fi
rm fem1d_heat_steady_prb
#
echo "Program output written to fem1d_heat_steady_prb_output.txt"
