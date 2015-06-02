#!/bin/bash
#
cp cube_felippa_rule.h /$HOME/include
#
gcc -c -I /$HOME/include cube_felippa_rule.c
if [ $? -ne 0 ]; then
  echo "Errors compiling cube_felippa_rule.c"
  exit
fi
#
mv cube_felippa_rule.o ~/libc/$ARCH/cube_felippa_rule.o
#
echo "Library installed as ~/libc/$ARCH/cube_felippa_rule.o"
