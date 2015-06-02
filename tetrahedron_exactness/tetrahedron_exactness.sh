#!/bin/bash
#
gcc -c tetrahedron_exactness.c
if [ $? -ne 0 ]; then
  echo "Errors compiling tetrahedron_exactness.c"
  exit
fi
#
gcc tetrahedron_exactness.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tetrahedron_exactness.o"
  exit
fi
rm tetrahedron_exactness.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/tetrahedron_exactness
#
echo "Program installed as ~/binc/$ARCH/tetrahedron_exactness"
