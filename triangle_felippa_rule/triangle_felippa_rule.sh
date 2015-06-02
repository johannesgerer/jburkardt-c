#!/bin/bash
#
cp triangle_felippa_rule.h /$HOME/include
#
gcc -c -I /$HOME/include triangle_felippa_rule.c
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_felippa_rule.c"
  exit
fi
#
mv triangle_felippa_rule.o ~/libc/$ARCH/triangle_felippa_rule.o
#
echo "Library installed as ~/libc/$ARCH/triangle_felippa_rule.o"
