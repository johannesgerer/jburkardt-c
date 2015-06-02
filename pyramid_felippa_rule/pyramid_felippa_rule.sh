#!/bin/bash
#
cp pyramid_felippa_rule.h /$HOME/include
#
gcc -c -I /$HOME/include pyramid_felippa_rule.c
if [ $? -ne 0 ]; then
  echo "Errors compiling pyramid_felippa_rule.c"
  exit
fi
#
mv pyramid_felippa_rule.o ~/libc/$ARCH/pyramid_felippa_rule.o
#
echo "Library installed as ~/libc/$ARCH/pyramid_felippa_rule.o"
