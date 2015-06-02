#!/bin/bash
#
cp square_symq_rule.h /$HOME/include
#
gcc -c -I/$HOME/include square_symq_rule.c
if [ $? -ne 0 ]; then
  echo "Errors compiling square_symq_rule.c"
  exit
fi
#
mv square_symq_rule.o ~/libc/$ARCH/square_symq_rule.o
#
echo "Library installed as ~/libc/$ARCH/square_symq_rule.o"
