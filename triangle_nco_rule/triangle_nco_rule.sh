#!/bin/bash
#
cp triangle_nco_rule.h /$HOME/include
#
gcc -c -I /$HOME/include triangle_nco_rule.c
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_nco_rule.c"
  exit
fi
#
mv triangle_nco_rule.o ~/libc/$ARCH/triangle_nco_rule.o
#
echo "Library installed as ~/libc/$ARCH/triangle_nco_rule.o"
