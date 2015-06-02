#!/bin/bash
#
gcc -c -g -I/$HOME/include fd2d_heat_steady_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fd2d_heat_steady_prb.c"
  exit
fi
rm compiler.txt
#
gcc fd2d_heat_steady_prb.o /$HOME/libc/$ARCH/fd2d_heat_steady.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fd2d_heat_steady_prb.o."
  exit
fi
#
rm fd2d_heat_steady_prb.o
#
mv a.out fd2d_heat_steady_prb
./fd2d_heat_steady_prb > fd2d_heat_steady_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fd2d_heat_steady_prb."
  exit
fi
rm fd2d_heat_steady_prb
#
echo "Program output written to fd2d_heat_steady_prb_output.txt"
