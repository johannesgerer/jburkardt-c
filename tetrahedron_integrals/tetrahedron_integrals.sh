#!/bin/bash
#
cp tetrahedron_integrals.h /$HOME/include
#
gcc -c -g -I /$HOME/include tetrahedron_integrals.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tetrahedron_integrals.c"
  exit
fi
rm compiler.txt
#
mv tetrahedron_integrals.o ~/libc/$ARCH/tetrahedron_integrals.o
#
echo "Library installed as ~/libc/$ARCH/tetrahedron_integrals.o"
