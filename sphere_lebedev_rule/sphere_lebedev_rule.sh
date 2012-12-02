#!/bin/bash
#
cp sphere_lebedev_rule.h /$HOME/include
#
gcc -c -g sphere_lebedev_rule.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_lebedev_rule.c."
  exit
fi
rm compiler.txt
#
mv sphere_lebedev_rule.o ~/libc/$ARCH/sphere_lebedev_rule.o
#
echo "Library installed as ~/libc/$ARCH/sphere_lebedev_rule.o"
