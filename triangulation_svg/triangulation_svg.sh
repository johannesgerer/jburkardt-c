#!/bin/bash
#
gcc -c triangulation_svg.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_svg.c."
  exit
fi
rm compiler.txt
#
gcc triangulation_svg.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_svg.o."
  exit
fi
#
rm triangulation_svg.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/triangulation_svg
#
echo "Executable installed as ~/binc/$ARCH/triangulation_svg"
