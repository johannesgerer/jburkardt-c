#!/bin/bash
#
gcc -c -I/$HOME/include triangle_svg_prb.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_svg_prb.c"
  exit
fi
rm compiler.txt
#
gcc triangle_svg_prb.o /$HOME/libc/$ARCH/triangle_svg.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangle_svg_prb.o."
  exit
fi
#
rm triangle_svg_prb.o
#
mv a.out triangle_svg_prb
./triangle_svg_prb > triangle_svg_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running triangle_svg_prb."
  exit
fi
rm triangle_svg_prb
#
echo "Program output written to triangle_svg_prb_output.txt"
