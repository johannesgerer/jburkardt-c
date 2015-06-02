#!/bin/bash
#
gcc -c tetrahedron_properties.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tetrahedron_properties.c"
  exit
fi
rm compiler.txt
#
g++ tetrahedron_properties.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tetrahedron_properties.o"
  exit
fi
rm tetrahedron_properties.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/tetrahedron_properties
#
echo "Executable installed as ~/binc/$ARCH/tetrahedron_properties"
