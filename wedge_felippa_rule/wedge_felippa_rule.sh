#!/bin/bash
#
cp wedge_felippa_rule.h /$HOME/include
#
gcc -c -I/$HOME/include wedge_felippa_rule.c
if [ $? -ne 0 ]; then
  echo "Errors compiling wedge_felippa_rule.c"
  exit
fi
#
mv wedge_felippa_rule.o ~/libc/$ARCH/wedge_felippa_rule.o
#
echo "Library installed as ~/libc/$ARCH/wedge_felippa_rule.o"
