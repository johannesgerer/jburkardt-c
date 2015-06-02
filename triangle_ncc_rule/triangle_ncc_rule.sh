#!/bin/bash
#
cp triangle_ncc_rule.h /$HOME/include
#
gcc -c -I /$HOME/include triangle_ncc_rule.c
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_ncc_rule.c"
  exit
fi
#
mv triangle_ncc_rule.o ~/libc/$ARCH/triangle_ncc_rule.o
#
echo "Library installed as ~/libc/$ARCH/triangle_ncc_rule.o"
