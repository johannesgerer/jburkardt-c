#!/bin/bash
#
cp tetrahedron_arbq_rule.h /$HOME/include
#
gcc -c -I/$HOME/include tetrahedron_arbq_rule.c
if [ $? -ne 0 ]; then
  echo "Errors compiling tetrahedron_arbq_rule.c"
  exit
fi
#
mv tetrahedron_arbq_rule.o ~/libc/$ARCH/tetrahedron_arbq_rule.o
#
echo "Library installed as ~/libc/$ARCH/tetrahedron_arbq_rule.o"
