#!/bin/bash
#
cp simplex_gm_rule.h /$HOME/include
#
gcc -c -I /$HOME/include simplex_gm_rule.c
if [ $? -ne 0 ]; then
  echo "Errors compiling simplex_gm_rule.c"
  exit
fi
#
mv simplex_gm_rule.o ~/libc/$ARCH/simplex_gm_rule.o
#
echo "Library installed as ~/libc/$ARCH/simplex_gm_rule.o"
