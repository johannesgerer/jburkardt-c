#!/bin/bash
#
cp square_arbq_rule.h /$HOME/include
#
gcc -c -I/$HOME/include square_arbq_rule.c
if [ $? -ne 0 ]; then
  echo "Errors compiling square_arbq_rule.c"
  exit
fi
#
mv square_arbq_rule.o ~/libc/$ARCH/square_arbq_rule.o
#
echo "Library installed as ~/libc/$ARCH/square_arbq_rule.o"
