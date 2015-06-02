#!/bin/bash
#
cp cube_arbq_rule.h /$HOME/include
#
gcc -c -I/$HOME/include cube_arbq_rule.c
if [ $? -ne 0 ]; then
  echo "Errors compiling cube_arbq_rule.c"
  exit
fi
#
mv cube_arbq_rule.o ~/libc/$ARCH/cube_arbq_rule.o
#
echo "Library installed as ~/libc/$ARCH/cube_arbq_rule.o"
