#!/bin/bash
#
gcc -c fem_to_xml.c
if [ $? -ne 0 ]; then
  echo "Errors compiling fem_to_xml.c"
  exit
fi
#
gcc fem_to_xml.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fem_to_xml.o."
  exit
fi
#
rm fem_to_xml.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/fem_to_xml
#
echo "Executable installed as ~/binc/$ARCH/fem_to_xml"
