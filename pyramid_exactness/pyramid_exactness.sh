#!/bin/bash
#
gcc -c pyramid_exactness.c
if [ $? -ne 0 ]; then
  echo "Errors compiling pyramid_exactness.c"
  exit
fi
#
gcc pyramid_exactness.o -lm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading pyramid_exactness.o"
  exit
fi
rm pyramid_exactness.o
#
chmod ugo+x a.out
mv a.out ~/binc/$ARCH/pyramid_exactness
#
echo "Executable installed as ~/binc/$ARCH/pyramid_exactness"
