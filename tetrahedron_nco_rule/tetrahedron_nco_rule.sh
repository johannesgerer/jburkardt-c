#!/bin/bash
#
cp tetrahedron_nco_rule.h /$HOME/include
#
gcc -c -I /$HOME/include tetrahedron_nco_rule.c
if [ $? -ne 0 ]; then
  echo "Errors compiling tetrahedron_nco_rule.c"
  exit
fi
#
mv tetrahedron_nco_rule.o ~/libc/$ARCH/tetrahedron_nco_rule.o
#
echo "Library installed as ~/libc/$ARCH/tetrahedron_nco_rule.o"
