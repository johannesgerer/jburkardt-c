#!/bin/bash
#
cp tetrahedron_ncc_rule.h /$HOME/include
#
gcc -c -I /$HOME/include tetrahedron_ncc_rule.c
if [ $? -ne 0 ]; then
  echo "Errors compiling tetrahedron_ncc_rule.c"
  exit
fi
#
mv tetrahedron_ncc_rule.o ~/libc/$ARCH/tetrahedron_ncc_rule.o
#
echo "Library installed as ~/libc/$ARCH/tetrahedron_ncc_rule.o"
