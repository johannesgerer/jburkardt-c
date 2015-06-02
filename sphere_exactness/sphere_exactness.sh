#!/bin/bash
#
gcc -c sphere_exactness.c
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_exactness.cpp"
  exit
fi
#
gcc sphere_exactness.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sphere_exactness.o"
  exit
fi
rm sphere_exactness.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/sphere_exactness
#
echo "Program installed as ~/binc/$ARCH/sphere_exactness"
