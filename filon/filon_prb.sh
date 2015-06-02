#!/bin/bash
#
gcc -c -I/$HOME/include filon_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling filon_prb.c"
  exit
fi
rm compiler.txt
#
gcc filon_prb.o /$HOME/libc/$ARCH/filon.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading filon_prb.o."
  exit
fi
#
rm filon_prb.o
#
mv a.out filon_prb
./filon_prb > filon_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running filon_prb."
  exit
fi
rm filon_prb
#
echo "Program output written to filon_prb_output.txt"
