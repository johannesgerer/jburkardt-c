#!/bin/bash
#
cp triangle_symq_rule.h /$HOME/include
#
gcc -c -I/$HOME/include triangle_symq_rule.c
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_symq_rule.c"
  exit
fi
#
mv triangle_symq_rule.o ~/libc/$ARCH/triangle_symq_rule.o
#
echo "Library installed as ~/libc/$ARCH/triangle_symq_rule.o"
