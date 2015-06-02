#!/bin/bash
#
gcc -c triangle_to_xml.c
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_to_xml.c"
  exit
fi
#
gcc triangle_to_xml.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangle_to_xml.o."
  exit
fi
#
rm triangle_to_xml.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/triangle_to_xml
#
echo "Executable installed as ~/binc/$ARCH/triangle_to_xml"
